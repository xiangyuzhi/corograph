#include "galois/AtomicHelpers.h"
#include "galois/Bag.h"
#include "galois/graphs/LCGraph.h"
#include "galois/runtime/Executor_EdgeMap.h"
#include "galois/substrate/ThreadPool.h"

unsigned int stepShift = 13;
std::string inputFile;
unsigned int reportNode = 4819611;
int numThreads = 1;
const uint32 PSIZE = 18;

using Graph = galois::graphs::graph<uint32>;
struct UpdateRequestIndexer {
  typedef std::less<uint32> compare;
  unsigned shift;
  template <typename R> unsigned int operator()(const R &req) const {
    if (req.dist < 10)
      return req.dist >> shift;
    return 10;
  }
};
constexpr static const unsigned CHUNK_SIZE = 1024U;
constexpr static const unsigned CG_CHUNK_SIZE = 4096U;

using vw = galois::graphs::vertex_warp<uint32>;
using pw = galois::graphs::part_wrap<uint32>;
namespace gwl = galois::worklists;
using PSchunk = gwl::CM<CHUNK_SIZE, vw>;
using SGchunk = gwl::CM2<CG_CHUNK_SIZE, pw>;
using Ck = gwl::CK<CHUNK_SIZE, vw>;
using Ck2 = gwl::CK<CG_CHUNK_SIZE, pw>;
using OBIM = gwl::OBIM<UpdateRequestIndexer, PSchunk, SGchunk, Ck, Ck2>;
typedef std::pair<uint32, uint32> label_type;

struct LNode {
  std::atomic<unsigned int> comp_current;
  unsigned int comp_old;
};

struct LNode2 {
  unsigned int comp_current;
  unsigned int comp_old;
};

struct CC_F {
  LNode2 *vdata; // distance
  galois::GReduceLogicalOr &change;
  explicit CC_F(LNode2 *_label, galois::GReduceLogicalOr &_c)
      : vdata(_label), change(_c) {}

  inline bool filterFunc(uint32 src) const {
    if (vdata[src].comp_old > vdata[src].comp_current) {
      vdata[src].comp_old = vdata[src].comp_current;
      change.update(true);
      return false;
    }
    return true;
  }
  inline uint32 scatterFunc(uint32 src) const {
    return vdata[src].comp_current;
  }
  inline bool gatherFunc(uint32 updateVal, uint32 destId) const {
    vdata[destId].comp_current =
        std::min(updateVal, vdata[destId].comp_current);
    return false;
  }
  inline uint32 pushFunc(uint32 dst, uint32 newVal) const { return dst; }
  static inline unsigned int applyWeight(unsigned int weight,
                                         unsigned int updateVal) {
    return updateVal;
  }
};

void sync_cc(Graph &graph, auto &tt, LNode2 *label) {
  galois::InsertBag<uint32> tmp;
  galois::GReduceLogicalOr changed;
  galois::runtime::syncExecutor exec(graph, CC_F(label, changed));
  uint32 iter = 0;
  do {
    printf("iter %d\n", ++iter);
    changed.reset();
    exec.EdgeMap(tt, tmp);
  } while (changed.reduce());
}

struct CC_F2 {
  uint32 *vdata; // distance
  explicit CC_F2(uint32 *_label) : vdata(_label) {}

  inline bool filterFunc(uint32 src, uint32 lb) const {
    return (vdata[src] < lb);
  }
  inline bool gatherFunc(uint32 updateVal, uint32 destId) const {
    if (updateVal < vdata[destId]) {
      vdata[destId] = updateVal;
      return true;
    }
    return false;
  }
  inline vw pushFunc(uint32 dst, uint32 newVal) const {
    return vw(dst, newVal);
  }
  static inline unsigned int applyWeight(unsigned int weight,
                                         unsigned int updateVal) {
    return updateVal;
  }
};

template <typename OBIMTy = OBIM>
void async_cc(Graph &graph, auto &initFrontier, uint32 *label) {
  galois::runtime::asyncPriorityEdgeMap<OBIMTy>(
      graph, UpdateRequestIndexer{stepShift}, CC_F2(label),
      galois::iterate(initFrontier));
}

void init_galois(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    if (i + 1 != argc) {
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
    printf("Usage : %s <filename> -s <start node(not needed for pr)> -t "
           "<numThreads(optional)> -iter <#iterations(optional) -rounds "
           "<#rounds(default 3)> \n",
           argv[0]);
    exit(1);
  }
  inputFile = std::string(argv[1]);
  numThreads = galois::setActiveThreads(numThreads);
}

template <typename TMP> struct Range {
  TMP &tmp;
  Range(TMP &t) : tmp(t) {}
  void operator()(unsigned tid, unsigned total) { tmp.range(tid, total); }
};
template <typename TMP>
void initRange(TMP &tmp) { // the second dispatch goto here
  Range<TMP> ranger(tmp);
  galois::on_each(ranger);
}

class temp : private galois::graphs::internal::LocalIteratorFeature<true> {
  uint32 num;

public:
  temp(uint32 _num) : num(_num) {}
  void range(uint32 tid, uint32 total) {
    uint32 len = num / total + 1;
    this->setLocalRange(len * tid, std::min(num, len * (tid + 1)));
  }
  using iterator = boost::counting_iterator<uint32>;
  typedef iterator local_iterator;
  local_iterator local_begin() { return iterator(this->localBegin(num)); }
  local_iterator local_end() { return iterator(this->localEnd(num)); }
  iterator begin() const { return iterator(0); }
  iterator end() const { return iterator(num); }
};

int main(int argc, char **argv) {
  galois::substrate::ThreadPool tp;
  std::unique_ptr<galois::substrate::internal::LocalTerminationDetection<>>
      m_termPtr;
  std::unique_ptr<galois::substrate::internal::BarrierInstance<>> m_biPtr;
  galois::substrate::internal::setThreadPool(
      &tp); // init the thread and set the pointer ----
  m_biPtr = std::make_unique<galois::substrate::internal::BarrierInstance<>>();
  m_termPtr = std::make_unique<
      galois::substrate::internal::LocalTerminationDetection<>>();
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
  galois::preAlloc(
      numThreads +
      approxNodeData /
          galois::runtime::pagePoolSize()); // preAlloc(nodesize/pagesize)

  auto *label = new uint32[G.numV];
  //    auto* label = new label_type [G.numV];
  //    auto* label = new LNode [G.numV];

  temp tt(G.numV);
  initRange(tt);

  for (uint32 _ = 0; _ < 5; _++) {
    galois::InsertBag<vw> initFrontier;
    galois::do_all(
        galois::iterate(tt),
        [&](const uint32 &n) {
          label[n] = n;
          initFrontier.push_back(vw(n, n));
          //                label[n].second = label[n].first = n;
        },
        galois::no_stats(), galois::loopname("initNodeData"));

    std::cout << "Running connected components algorithm\n";

    struct timespec start, end;
    float time;
    if (clock_gettime(CLOCK_REALTIME, &start) == -1) {
      perror("clock gettime");
    }
    //        sync_cc(G, tt, label);
    async_cc(G, initFrontier, label);
    if (clock_gettime(CLOCK_REALTIME, &end) == -1) {
      perror("clock gettime");
    }
    time =
        (end.tv_sec - start.tv_sec) + (int)(end.tv_nsec - start.tv_nsec) / 1e9;
    printf("time: %lf sec\n", time);

    //        std::cout << "Node " << reportNode << " has label " <<
    //        label[report].comp_current << "\n";
    std::cout << "Node " << reportNode << " has label " << label[report]
              << "\n";
    uint32_t ccnt = 0;
    std::map<uint32, uint32> mm;
    for (uint32_t i = 0; i < G.numV; i++) {
      //            mm[label[i].comp_current] += 1;
      mm[label[i]] += 1;
    }
    printf("component num: %d \n", mm.size());
  }

  return 0;
}
