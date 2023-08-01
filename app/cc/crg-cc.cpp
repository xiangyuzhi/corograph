#include "galois/substrate/ThreadPool.h"
#include "galois/graphs/LCGraph.h"
#include "galois/Bag.h"
#include "galois/runtime/Executor_EdgeMap.h"
#include "galois/AtomicHelpers.h"

unsigned int stepShift = 13;
std::string inputFile;
unsigned int reportNode = 4819611;
int numThreads = 1;
const uint32 PSIZE = 18;

using Graph = galois::graphs::graph<uint32>;
struct UpdateRequestIndexer {
    typedef std::less<uint32> compare;
    unsigned shift;
    template <typename R>
    unsigned int operator()(const R& req) const {
        if(req.dist < 10) return req.dist>>shift;
        return 10;
    }
};
constexpr static const unsigned CHUNK_SIZE      = 1024U;
constexpr static const unsigned CG_CHUNK_SIZE   = 4096U;


using vw = galois::graphs::vertex_warp<uint32>;
using pw = galois::graphs::part_wrap<uint32>;
namespace gwl = galois::worklists;
using PSchunk = gwl::CM<CHUNK_SIZE, vw>;
using SGchunk = gwl::CM2<CG_CHUNK_SIZE, pw>;
using Ck      = gwl::CK<CHUNK_SIZE, vw>;
using Ck2     = gwl::CK<CG_CHUNK_SIZE, pw>;
using OBIM    = gwl::OBIM<UpdateRequestIndexer, PSchunk, SGchunk, Ck, Ck2>;
typedef std::pair<uint32, uint32> label_type;

struct LNode {
    std::atomic<unsigned int> comp_current;
    unsigned int comp_old;
};

struct LNode2 {
    unsigned int comp_current;
    unsigned int comp_old;
};


void cc(Graph& graph, auto& tt, LNode* label) {
    galois::GReduceLogicalOr changed;
    uint32 iter = 0;
    do {
        printf("iter %d\n", ++iter);
        changed.reset();
        galois::do_all(
            galois::iterate(tt),
            [&](const uint32& src) {
                LNode& sdata = label[src];//graph.getData(src, galois::MethodFlag::UNPROTECTED);
                if (sdata.comp_old > sdata.comp_current) {
                    sdata.comp_old = sdata.comp_current;
                    unsigned int label_new = sdata.comp_current;
                    changed.update(true);
                    for(uint32 e =graph.offset[src]; e<graph.offset[src+1];e++) {
                        uint32 dst             = graph.ngh[e];
                        auto& ddata            = label[dst];
                        galois::atomicMin(ddata.comp_current, label_new);
                    }
                }
            },
            galois::disable_conflict_detection(), galois::steal(),
            galois::loopname("LabelPropAlgo"));
    } while (changed.reduce());
}

void init_galois(int argc, char** argv) {
    for (int i = 1; i < argc; i++){
        if (i + 1 != argc){
            if (strcmp(argv[i], "-delta") == 0) { // delta
                stepShift = (unsigned int)atoi(argv[i + 1]);
                i++;
            }
            if (strcmp(argv[i], "-t") == 0) { // num threads
                numThreads = (int)atoi(argv[i + 1]);
                i++;
            }
        }
    }
    if (argc < 2) {
        printf("Usage : %s <filename> -s <start node(not needed for pr)> -t <numThreads(optional)> -iter <#iterations(optional) -rounds <#rounds(default 3)> \n", argv[0]);
        exit(1);
    }
    inputFile = std::string(argv[1]);
    numThreads = galois::setActiveThreads(numThreads);
}


template <typename TMP>
struct Range {
    TMP& tmp;
    Range(TMP& t) : tmp(t){}
    void operator()(unsigned tid, unsigned total) {
        tmp.range(tid, total);
    }
};
template <typename TMP>
void readGraphDispatch(TMP& tmp) { // the second dispatch goto here
    Range<TMP> ranger(tmp);
    galois::on_each(ranger);
}

class temp : private galois::graphs::internal::LocalIteratorFeature<true>{
    uint32 num;
public:
    temp(uint32 _num): num(_num) {}
    void range(uint32 tid, uint32 total){
        uint32 len = num/total+1;
        this->setLocalRange(len*tid, std::min(num,len*(tid+1)) );
    }
    using iterator = boost::counting_iterator<uint32>;
    typedef iterator local_iterator;
    local_iterator local_begin() {
        return iterator(this->localBegin(num));
    }
    local_iterator local_end() {
        return iterator(this->localEnd(num));
    }
    iterator begin() const { return iterator(0); }
    iterator end() const { return iterator(num); }
};

int main(int argc, char** argv) {
    galois::substrate::ThreadPool tp;
    std::unique_ptr<galois::substrate::internal::LocalTerminationDetection<>> m_termPtr;
    std::unique_ptr<galois::substrate::internal::BarrierInstance<>> m_biPtr;
    galois::substrate::internal::setThreadPool(&tp); // init the thread and set the pointer ----
    m_biPtr   = std::make_unique<galois::substrate::internal::BarrierInstance<>>();
    m_termPtr = std::make_unique<galois::substrate::internal::LocalTerminationDetection<>>();
    galois::substrate::internal::setBarrierInstance(m_biPtr.get());
    galois::substrate::internal::setTermDetect(m_termPtr.get());
    galois::runtime::internal::PageAllocState<> m_pa;
    galois::runtime::internal::setPagePoolState(&m_pa);
    init_galois(argc, argv);
    Graph G;
    commandLine P(argc, argv);

    uint32 report = reportNode;

    galois::graphs::init_graph(G, P);

    std::cout << "Read " << G.numV << " nodes, " << G.numE << " edges\n";

    printf("Partition Graph\n");
    partition(G, PSIZE);

    size_t approxNodeData = G.numV * 64;
    galois::preAlloc(numThreads + approxNodeData / galois::runtime::pagePoolSize());// preAlloc(nodesize/pagesize)

    auto* label = new LNode [G.numV];

    temp tt(G.numV);
    readGraphDispatch(tt);

    for(uint32 _=0; _<5;_++){
        galois::InsertBag<vw> initFrontier;
        galois::do_all(
            galois::iterate(tt),
            [&](const uint32& n) {
                label[n].comp_current = n;label[n].comp_old = MAX_NUM;
            },
            galois::no_stats(), galois::loopname("initNodeData"));

        std::cout << "Running connected components algorithm\n";

        struct timespec start, end;
        float time; if( clock_gettime( CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
        cc(G, tt, label);
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");} time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("time: %lf sec\n",time);

        std::cout << "Node " << reportNode << " has label " << label[report].comp_current << "\n";
        uint32_t ccnt =0;
        std::map<uint32, uint32> mm;
        for(uint32_t i=0;i<G.numV; i++) {
            mm[label[i].comp_current] += 1;
        }
        printf("component num: %d \n", mm.size());
    }


    return 0;
}


