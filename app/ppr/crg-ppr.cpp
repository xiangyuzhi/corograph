#include "galois/substrate/ThreadPool.h"
#include "galois/AtomicHelpers.h"
#include "galois/graphs/LCGraph.h"

std::string inputFile;
unsigned int startNode = 9;
unsigned int reportNode = 9;
int numThreads = 1;
const uint32 PSIZE = 18;

using Graph = galois::graphs::graph<uint32>;
struct UpdateRequestIndexer {
    typedef std::greater<uint32> compare;
    template <typename R>
    unsigned int operator()(const R& req) const {
//        unsigned int t = req.dist*1e9;
//        return t/1e6;
        return 1;
    }
};
constexpr static const unsigned CHUNK_SIZE      = 1024U;
constexpr static const unsigned CG_CHUNK_SIZE   = 4096U;


using vw = galois::graphs::vertex_warp<float>;
using pw = galois::graphs::part_wrap<float>;
using pd = std::pair<float, float>;
namespace gwl = galois::worklists;
using PSchunk = gwl::CM<CHUNK_SIZE, vw>;
using SGchunk = gwl::CM2<CG_CHUNK_SIZE, pw>;
using Ck      = gwl::CK<CHUNK_SIZE, vw>;
using Ck2     = gwl::CK<CG_CHUNK_SIZE, pw>;
using OBIM    = gwl::OBIM<UpdateRequestIndexer, PSchunk, SGchunk, Ck, Ck2>;


const float alpha = 0.01;
const float epsilon = 0.0000001;
const float k1 = 2 * alpha / (1 + alpha);
const float k2 = (1 - alpha) / (1 + alpha);

struct PPR_F{
    pd* vdata; // distance
    uint32* deg;
    explicit PPR_F(pd* _pr, uint32* _deg) : vdata(_pr), deg(_deg) {}

    inline bool filterFunc(uint32 src, float val) const{
        return false;
    }
    inline bool gatherFunc (float updateVal, uint32 destId) const{
        float oldval = vdata[destId].second;
        float epsd = epsilon * deg[destId];
        vdata[destId].second += updateVal;
        if(oldval < epsd && epsd < vdata[destId].second){
            return true;
        }
        return false;
    }
    inline vw pushFunc(uint32 dst, float newpr) const {
        vdata[dst].first += k1 * vdata[dst].second;
        float val = vdata[dst].second * k2/ deg[dst];
        vdata[dst].second = 0;
        return vw(dst, val);
    }
    static inline float applyWeight (unsigned int weight, float updateVal){
        return updateVal;
    }
};


template <typename OBIMTy = OBIM>
void parallel_nibble(Graph& graph, uint32 source, pd* pr) {

    PPR_F pf(pr, graph.deg);
    std::deque<vw> initFrontier;
    initFrontier.push_back(pf.pushFunc(source, 1.0));

    galois::runtime::asyncPriorityEdgeMap<OBIMTy>(
            graph,
            initFrontier,
            UpdateRequestIndexer{},
            pf);
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
    partition(G, PSIZE);

    size_t approxNodeData = G.numV * 64;
    galois::preAlloc(numThreads + approxNodeData / galois::runtime::pagePoolSize());// preAlloc(nodesize/pagesize)

    auto* pr = new pd [G.numV];

    for(uint32 _=0; _<5;_++){
        for(uint32 i=0;i<G.numV;i++){
            pr[i] = std::make_pair(0.0, 0.0);
        }
        pr[source] = std::make_pair(0.0, 1.0);

        std::cout << "Running parallel nibble algorithm\n";
        struct timespec start, end;float time; if( clock_gettime( CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
        parallel_nibble(G, source, pr);
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");} time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;printf("time: %lf sec\n",time);
        std::cout << "Node " << reportNode << " has ppr " << pr[report].first << "\n";

        float maxpr =0.0;
        for(uint32_t i=0;i<G.numV; i++) {
            maxpr = std::max(maxpr, pr[i].first);
        }
        printf("max ppr: %.8f \n",maxpr);
    }


    return 0;
}


