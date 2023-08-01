#include "galois/substrate/ThreadPool.h"
#include "galois/AtomicHelpers.h"
#include "galois/graphs/LCGraph.h"
#include "galois/Bag.h"

//unsigned int corenum = 4;
std::string inputFile;
unsigned int startNode = 9;
unsigned int reportNode = 4819611;
int numThreads = 1;
const uint32 PSIZE = 18;

using Graph = galois::graphs::graph<uint32>;
struct UpdateRequestIndexer {
    typedef std::greater<uint32> compare;
    template <typename R>
    unsigned int operator()(const R& req) const {
        return 1;
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

struct KCORE_F{
    unsigned int* vdata; // curdeg
    uint32 corenum;
    explicit KCORE_F(unsigned int* _curdeg, uint32 _cn) : vdata(_curdeg), corenum(_cn){}

    inline bool filterFunc(uint32 src, uint32 dis) const{
        return false;
    }
    inline bool gatherFunc (unsigned int updateVal, uint32 destId) const{
        vdata[destId] --;
        if(vdata[destId] == corenum - 1){
            return true;
        }
        return false;
    }
    inline vw pushFunc(uint32 dst, uint32 newdis) const {
        return vw(dst, vdata[dst]);
    }
    static inline unsigned int applyWeight (unsigned int updateVal, unsigned int weight){
        return 0;
    }
};


template <typename OBIMTy = OBIM>
void kcore(Graph& graph, auto &initFrontier, uint32* curdeg, uint32 corenum) {

    galois::runtime::asyncPriorityEdgeMap<OBIMTy>(
            graph,
            UpdateRequestIndexer{},
            KCORE_F(curdeg, corenum),
            galois::iterate(initFrontier));
}


void init_galois(int argc, char** argv) {
    for (int i = 1; i < argc; i++){
        if (i + 1 != argc){
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

    uint32 source = startNode;
    uint32 report = reportNode;

    galois::graphs::init_graph(G, P);

    std::cout << "Read " << G.numV << " nodes, " << G.numE << " edges\n";

    printf("Partition Graph\n");
    partition(G, numThreads);
    temp tt(G.numV);
    readGraphDispatch(tt);

    size_t approxNodeData = G.numV * 64;
    galois::preAlloc(numThreads + approxNodeData / galois::runtime::pagePoolSize());// preAlloc(nodesize/pagesize)

    auto* curdeg = new uint32 [G.numV];
    auto* ans = new uint32 [G.numV];


    for(uint32 _=0; _<1;_++){

        std::cout << "Running sorted k-core algorithm\n";
        galois::do_all(
                galois::iterate(tt),
                [&](const uint32& n) {
                    curdeg[n] = G.deg[n];
                },galois::no_stats(), galois::loopname("initNodeData"));
        galois::GReduceLogicalOr changed;
        struct timespec start, end;
        float time; if( clock_gettime( CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
        for(uint32 corenum =1;corenum <G.numV;corenum++){
            galois::InsertBag<vw> initFrontier;
            galois::do_all(
                    galois::iterate(tt),
                    [&](const uint32& n) {
                        if(curdeg[n] == corenum-1){
                            initFrontier.emplace_back(n, curdeg[n]);
                        }
                    },
                    galois::no_stats(), galois::loopname("initNodeData"));
            kcore(G, initFrontier, curdeg, corenum);
            changed.reset();
            galois::do_all(
                    galois::iterate(tt),
                    [&](const uint32& n) {
                        if(curdeg[n] >= corenum ){
                            ans[n] = corenum;
                            changed.update(true);
                        }
                    },galois::no_stats(), galois::loopname("initNodeData"));
            if(!changed.reduce()) {
                printf("max core number: %d\n",corenum-1);
                break;
            }
        }
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");} time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("time: %lf sec\n",time);
    }


    return 0;
}


