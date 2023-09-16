#include "galois/AtomicHelpers.h"
#include "galois/Bag.h"
#include "galois/graphs/LCGraph.h"
#include "galois/runtime/Executor_EdgeMap.h"
#include "galois/substrate/ThreadPool.h"

unsigned int stepShift = 13;
unsigned int startNode = 9;
unsigned int reportNode = 9;
int numThreads = 1;
const uint32 PSIZE = 18;

using Graph = galois::graphs::graph<uint32>;
struct UpdateRequestIndexer {
  typedef std::greater<uint32> compare;
  template <typename R> unsigned int operator()(const R &req) const {
    return 1;
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

const float alpha = 0.15;
const float epsilon = 0.000001;

struct PR_F {
  float *curpr;
  float *nextpr; // distance
  uint32 *deg;
  explicit PR_F(float *_cpr, float *_npr, uint32 *_deg)
      : curpr(_cpr), nextpr(_npr), deg(_deg) {}

  inline bool filterFunc(uint32 src) const { return false; }
  inline uint32 scatterFunc(uint32 src) const { return curpr[src] / deg[src]; }
  inline bool gatherFunc(float updateVal, uint32 destId) const {
    nextpr[destId] += updateVal;
    return false;
  }
  inline uint32 pushFunc(uint32 dst, float newpr) const { return dst; }
  static inline float applyWeight(unsigned int weight, float updateVal) {
    return updateVal;
  }
};

template <typename OBIMTy = OBIM> void pr(Graph &graph, auto &all, PR_F &pr) {

  galois::InsertBag<uint32> Frontier, nextF;
  galois::runtime::syncExecutor exec(graph, pr);
  exec.EdgeMap(all, nextF);
  galois::do_all(
      galois::iterate(all),
      [&](const uint32 &n) {
        pr.nextpr[n] = 0.15 / graph.numV + (0.85 * pr.nextpr[n]);
        if (std::abs(pr.nextpr[n] - pr.curpr[n]) > epsilon) {
          Frontier.push_back(n);
          pr.curpr[n] = 0.0;
        }
      },
      galois::no_stats(), galois::loopname("Reset"));
  for (uint32 _ = 0; _ < 9; _++) {
    exec.EdgeMap(Frontier, nextF);
    galois::do_all(
        galois::iterate(all),
        [&](const uint32 &n) {
          pr.nextpr[n] = 0.15 / graph.numV + (0.85 * pr.nextpr[n]);
          if (std::abs(pr.nextpr[n] - pr.curpr[n]) > epsilon) {
            Frontier.push_back(n);
            pr.curpr[n] = 0.0;
          }
        },
        galois::no_stats(), galois::loopname("Reset"));
    std::swap(pr.nextpr, pr.curpr);
  }
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

  uint32 source = startNode;
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

  temp all(G.numV);
  initRange(all);

  auto *curprv = new float[G.numV];
  auto *nextprv = new float[G.numV];
  PR_F prf(curprv, nextprv, G.deg);

  std::cout << "INFO: Using delta-step of " << (1 << stepShift) << "\n";

  for (uint32 _ = 0; _ < 5; _++) {
    galois::do_all(
        galois::iterate(all),
        [&](const uint32 &n) { curprv[n] = 1.0 / G.numV; }, galois::no_stats(),
        galois::loopname("Reset"));

    std::cout << "Running parallel nibble algorithm\n";
    struct timespec start, end;
    float time;
    if (clock_gettime(CLOCK_REALTIME, &start) == -1) {
      perror("clock gettime");
    }
    pr(G, all, prf);
    if (clock_gettime(CLOCK_REALTIME, &end) == -1) {
      perror("clock gettime");
    }
    time =
        (end.tv_sec - start.tv_sec) + (int)(end.tv_nsec - start.tv_nsec) / 1e9;
    printf("time: %lf sec\n", time);
    std::cout << "Node " << reportNode << " has pr " << curprv[report] << "\n";

    float maxpr = 0.0;
    for (uint32_t i = 0; i < G.numV; i++) {
      maxpr = std::max(maxpr, curprv[i]);
    }
    printf("max pr: %.8f \n", maxpr);
  }

  return 0;
}
