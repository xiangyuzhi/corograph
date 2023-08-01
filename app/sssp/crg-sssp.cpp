#include "galois/substrate/ThreadPool.h"
#include "galois/graphs/LCGraph.h"
#include "galois/Bag.h"

unsigned int stepShift = 13;
std::string inputFile;
unsigned int startNode = 9;
unsigned int reportNode = 4819611;
int numThreads = 1;


using Graph = galois::graphs::graph<uint32>;
struct UpdateRequestIndexer {
    typedef std::less<uint32> compare;
    unsigned shift;
    template <typename R>
    unsigned int operator()(const R& req) const {
        unsigned int t = req.dist >> shift;
        return t;
    }
};
constexpr static const unsigned CHUNK_SIZE      = 512U;
constexpr static const unsigned CG_CHUNK_SIZE   = 1024U;


using vw = galois::graphs::vertex_warp<uint32>;
//using pw = std::pair<uint32, uint32>;
using pw = galois::graphs::part_wrap<uint32>;
namespace gwl = galois::worklists;
using PSchunk = gwl::CM<CHUNK_SIZE, vw>;
using SGchunk = gwl::CM2<CG_CHUNK_SIZE, pw>;
using Ck      = gwl::CK<CHUNK_SIZE, vw>;
using Ck2     = gwl::CK<CG_CHUNK_SIZE, pw>;
using OBIM    = gwl::OBIM<UpdateRequestIndexer, PSchunk, SGchunk, Ck, Ck2>;

struct SSSP_F{
    unsigned int* vdata; // distance
    explicit SSSP_F(unsigned int* _distance):vdata(_distance){}

    inline bool filterFunc(uint32 src, uint32 dis) const{
        return vdata[src] < dis;
    }
    inline bool gatherFunc (unsigned int updateVal, uint32 destId) const{
        if(updateVal < vdata[destId]){
            vdata[destId] = updateVal;
            return true;
        }
        return false;
    }
    inline vw pushFunc(uint32 dst, uint32 newdis) const {
        return vw(dst, newdis);
    }
    static inline unsigned int applyWeight (unsigned int weight, unsigned int updateVal){
        return updateVal + weight;
    }
};


template <typename OBIMTy = OBIM>
void deltaStepAlgo(Graph& graph, auto& initFrontier, uint32* dist) {
    galois::runtime::asyncPriorityEdgeMap<OBIMTy>(
            graph,
            UpdateRequestIndexer{stepShift},
            SSSP_F(dist),
            galois::iterate(initFrontier));
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

    size_t approxNodeData = G.numV * 256;
    galois::preAlloc(numThreads + approxNodeData / galois::runtime::pagePoolSize());// preAlloc(nodesize/pagesize)

    auto* distance = new uint32 [G.numV];

    std::cout << "INFO: Using delta-step of " << (1 << stepShift) << "\n";

    for(uint32 _=0; _<5;_++){
        for(uint32 i=0;i<G.numV;i++){
            distance[i] = MAX_NUM;
        }
        distance[startNode] = 0;
        galois::InsertBag<vw> initFrontier;
        initFrontier.push_back(vw(source, 0));

        std::cout << "Running delta-step SSSP algorithm\n";
        struct timespec start, end;
        float time; if( clock_gettime( CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
        deltaStepAlgo(G, initFrontier, distance);
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");} time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("time: %lf sec\n",time);
        std::cout << "Node " << reportNode << " has distance " << distance[report] << "\n";

        uint32_t maxdist =0;
        for(uint32_t i=0;i<G.numV; i++) {
            if(distance[i]!=MAX_NUM)
                maxdist = std::max(maxdist, distance[i]);
        }
        printf("max distance: %d \n",maxdist);
    }

    return 0;
}


