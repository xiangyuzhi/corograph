#ifndef COROGRAPH_EXECUTOR_FOREACH_H
#define COROGRAPH_EXECUTOR_FOREACH_H

#include <algorithm>
#include <functional>
#include <memory>
#include <utility>

#include "galois/runtime/Context.h"
#include "galois/runtime/Corobj.h"

namespace galois {
//! Internal Galois functionality - Use at your own risk.
namespace runtime {

template <class WorkListTy, class Func, class Indexer, class Graph>
class PriorityExecutor {
protected:
  typedef typename WorkListTy::value_type value_type;
  typedef typename WorkListTy::part_type pw;
  struct ThreadLocalData {
    std::deque<value_type> facing;
    std::vector<std::vector<pw>> facing2;
    std::vector<galois::optional<pw>> low, high;
    std::vector<galois::optional<value_type>> tmp;

    Corobj<bool> coro_sca, coro_gat, coro_syn;
    explicit ThreadLocalData(uint32_t n) : facing(), facing2(n) {
      for (uint32_t i = 0; i < n;
           i++) { // seem not useful, the bottle neck is still memory bound
        facing2[i].reserve(100000);
      }
    }
  };

  substrate::TerminationDetection &term;
  WorkListTy wl;
  Func func;
  Graph &graph;

  void doScatter(auto *&p, auto &tld, bool &didWork) {
    while (p) {
      didWork = true;
      galois::optional<value_type> item;
      for (int i = 0; i < p->count; i++) { // filter
        item = p->extract_at(i);
        if (func.filterFunc((*item).vid, (*item).dist)) {
          continue;
        }
        // wl.update_cnt("edgecnt", graph.deg[(*item).vid]);
        tld.tmp.emplace_back(item);
      }
      p->clear();
      while (!tld.coro_sca())
        tld.coro_syn();
      tld.tmp.clear();
      p = wl.pop3();
    }
  }

  void doSync(auto &tld) {
    wl.scatter(tld.facing2);
    wl.sync();
  }

  void doGather(auto &tld) {
    while (wl.pop_part()) {
      auto *pg = wl.pop_gather();
      galois::optional<pw> it;
      while (pg) {
        for (int i = 0; i < pg->count; i++) {
          it = pg->extract_at(i);
          if (!((*it).e >> MAX_OFS))
            tld.low.emplace_back(it);
          else
            tld.high.emplace_back(it);
        }
        pg->clear();
        while (!tld.coro_gat()) {
        };
        tld.low.clear();
        tld.high.clear();
        pg = wl.pop_gather();
      }
      if (!tld.facing.empty()) {
        wl.push(tld.facing.begin(), tld.facing.end());
        tld.facing.clear();
      }
    }
  }

  bool runQueueSimple2(ThreadLocalData &tld) { // perthread execute
    bool didWork = false;
    auto *p = wl.pop2(); // pop a chunk
    while (p) {
      auto start = std::chrono::high_resolution_clock::now();
      doScatter(p, tld, didWork);
      auto end1 = std::chrono::high_resolution_clock::now();
      doSync(tld);
      auto end2 = std::chrono::high_resolution_clock::now();
      doGather(tld);
      auto end3 = std::chrono::high_resolution_clock::now();
      p = wl.slowPop3(); // pop the next priority queue
      auto end4 = std::chrono::high_resolution_clock::now();
      wl.updateTime("1Scatter", start, end1);
      wl.updateTime("2Sync", end1, end2);
      wl.updateTime("3Gather", end2, end3);
      wl.updateTime("4Pop", end3, end4);
    }
    return didWork;
  }

  template <bool isLeader> void go() {
    ThreadLocalData tld(graph.pnum);
    tld.coro_gat = coro_gather(tld.facing, tld.low, tld.high);
    tld.coro_sca = coro_scatter(tld.facing2, tld.tmp);
    tld.coro_syn = coro_sync(tld.facing2);
    do {
      bool b = runQueueSimple2(tld);
      term.localTermination(b);
    } while (!term.globalTermination());
  }

  Corobj<bool> coro_scatter(auto &ctx2, auto &tmp) {
    co_yield false;
    for (;;) {
      for (int id = 0; id < (int)tmp.size(); id += 64) {
        for (int prid = id; prid < std::min(id + 64, (int)tmp.size()); prid++) {
          _mm_prefetch(&graph.plgraph[(*tmp[prid]).vid], _MM_HINT_T0);
        }
        co_yield false;
        for (int prid = id; prid < std::min(id + 64, (int)tmp.size());
             prid++) { // scatter
          auto &vtxA = graph.plgraph[(*tmp[prid]).vid];
          auto *Arr = vtxA.PE;
          if (vtxA.deg2)
            _mm_prefetch(&graph.pledge[vtxA.offset], _MM_HINT_T0);
          for (uint32 e = 0; e < vtxA.deg1; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              ctx2[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                              (*tmp[prid]).dist);
            else
              ctx2[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*tmp[prid]).dist);
          }
          Arr = graph.pledge;
          for (uint32 e = vtxA.offset; e < vtxA.offset + vtxA.deg2; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              ctx2[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                              (*tmp[prid]).dist);
            else
              ctx2[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*tmp[prid]).dist);
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coro_gather(auto &ctx, auto &low, auto &high) {
    co_yield false;
    for (;;) {
      for (int id = 0; id < (int)low.size(); id += 64) {
        for (int prid = id; prid < std::min(id + 64, (int)low.size()); prid++) {
          _mm_prefetch(&func.vdata[(*low[prid]).e], _MM_HINT_T0);
        }
        co_yield false;
        for (int prid = id; prid < std::min(id + 64, (int)low.size()); prid++) {
          uint32 dst = (*low[prid]).e;
          auto update = func.applyWeight((*low[prid]).w, (*low[prid]).val);
          if (func.gatherFunc(update, dst)) {
            ctx.push_back(func.pushFunc(dst, update));
          }
        }
      }
      for (int id = 0; id < (int)high.size(); id += 64) {
        for (int prid = id; prid < std::min(id + 64, (int)high.size());
             prid++) {
          _mm_prefetch(&graph.highedge[(*high[prid]).w], _MM_HINT_T0);
        }
        co_yield false;
        for (int prid = id; prid < std::min(id + 64, (int)high.size());
             prid++) {
          uint32 deg = (*high[prid]).e & DW_HALF;
          auto dis = (*high[prid]).val;
          for (uint32 e = (*high[prid]).w; e < (*high[prid]).w + deg; e += 2) {
            uint32 dst = graph.highedge[e]; // here
            auto update = func.applyWeight(graph.highedge[e + 1], dis);
            if (func.gatherFunc(update, dst)) {
              ctx.push_back(func.pushFunc(dst, update));
            }
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coro_sync(auto &ctx2) {
    co_yield false;
    for (;;) {
      wl.scatter(ctx2);
      co_yield true;
    }
  }

public:
  PriorityExecutor(Func _func, Indexer indexer, Graph &_graph)
      : term(substrate::getSystemTermination(activeThreads)),
        wl(_graph.pnum, indexer), func(_func), graph(_graph) {}

  template <typename Frontier> void initFrontier(Frontier &front) {
    for (auto &a : front) {
      wl.push(a);
    }
    front.clear();
  }

  void print_profiling() {
    wl.summary_cnt();
    wl.summaryTime();
  }

  template <typename RangeTy> void initThread(const RangeTy &range) {
    auto rp = range.local_pair();
    wl.push(rp.first, rp.second);
    term.initializeThread();
  }

  void initThread() { term.initializeThread(); }

  void operator()() {
    bool isLeader = substrate::ThreadPool::isLeader();
    if (isLeader)
      go<true>();
    else
      go<false>();
  }
};

template <typename OBIM, typename Graph, typename Indexer, typename Func,
          typename Range> // RangeTy = LocalRange<InsertBag>
void asyncPriorityEdgeMap(Graph &graph, const Indexer &indexer, Func func,
                          const Range &range) {
  auto &barrier = getBarrier(activeThreads);
  PriorityExecutor<OBIM, Func, Indexer, Graph> W(func, indexer, graph);
  substrate::getThreadPool().run(
      activeThreads, [&W, &range]() { W.initThread(range()); },
      std::ref(barrier), std::ref(W));
  //            W.print_profiling();
}

template <class WorkListTy, class Func, class Indexer, class Graph>
class Executor {
protected:
  typedef typename WorkListTy::value_type value_type;
  typedef typename WorkListTy::part_type pw;
  struct ThreadLocalData {
    std::deque<value_type> facing;
    std::vector<std::vector<pw>> facing2;
    std::vector<galois::optional<pw>> low, high;
    std::vector<galois::optional<value_type>> tmp;

    Corobj<bool> coro_sca, coro_gat, coro_syn;
    explicit ThreadLocalData(uint32_t n) : facing(), facing2(n) {
      for (uint32_t i = 0; i < n;
           i++) { // seem not useful, the bottle neck is still memory bound
        facing2[i].reserve(100000);
      }
    }
  };

  substrate::TerminationDetection &term;
  WorkListTy wl;
  Func func;
  Graph &graph;

  void doScatter(auto *&p, auto &tld, bool &didWork) {
    while (p) {
      didWork = true;
      galois::optional<value_type> item;
      for (int i = 0; i < p->count; i++) { // filter
        item = p->extract_at(i);
        if (func.filterFunc((*item).vid, (*item).dist)) {
          continue;
        }
        tld.tmp.emplace_back(item);
      }
      p->clear();
      while (!tld.coro_sca())
        tld.coro_syn();
      tld.tmp.clear();
      p = wl.pop();
    }
  }

  void doSync(auto &tld) {
    wl.scatter(tld.facing2);
    wl.sync();
  }

  void doGather(auto &tld) {
    while (wl.pop_part()) {
      auto *pg = wl.pop_gather();
      galois::optional<pw> it;
      while (pg) {
        for (int i = 0; i < pg->count; i++) {
          it = pg->extract_at(i);
          if (!((*it).e >> MAX_OFS))
            tld.low.emplace_back(it);
          else
            tld.high.emplace_back(it);
        }
        pg->clear();
        while (!tld.coro_gat()) {
          if (!tld.facing.empty()) {
            wl.push(tld.facing.begin(), tld.facing.end());
            tld.facing.clear();
          }
        };
        tld.low.clear();
        tld.high.clear();
        pg = wl.pop_gather();
      }
    }
  }

  bool runQueueSimple2(ThreadLocalData &tld) { // perthread execute
    bool didWork = false;
    auto *p = wl.pop2(); // pop a chunk
    while (p) {
      auto start = std::chrono::high_resolution_clock::now();
      doScatter(p, tld, didWork);
      auto end1 = std::chrono::high_resolution_clock::now();
      doSync(tld);
      auto end2 = std::chrono::high_resolution_clock::now();
      doGather(tld);
      auto end3 = std::chrono::high_resolution_clock::now();
      p = wl.sloPop3(); // pop the next priority queue
      wl.updateTime("1Scatter", start, end1);
      wl.updateTime("2Sync", end1, end2);
      wl.updateTime("3Gather", end2, end3);
    }
    return didWork;
  }

  template <bool isLeader> void go() {
    ThreadLocalData tld(graph.pnum);
    tld.coro_gat = coro_gather(tld.facing, tld.low, tld.high);
    tld.coro_sca = coro_scatter(tld.facing2, tld.tmp);
    tld.coro_syn = coro_sync(tld.facing2);
    do {
      bool b = runQueueSimple2(tld);
      term.localTermination(b);
    } while (!term.globalTermination());
  }

  Corobj<bool> coro_scatter(auto &ctx2, auto &tmp) {
    co_yield false;
    for (;;) {
      for (int id = 0; id < (int)tmp.size(); id += 64) {
        for (int prid = id; prid < std::min(id + 64, (int)tmp.size()); prid++) {
          _mm_prefetch(&graph.plgraph[(*tmp[prid]).vid], _MM_HINT_T0);
        }
        co_yield false;
        for (int prid = id; prid < std::min(id + 64, (int)tmp.size());
             prid++) { // scatter
          auto &vtxA = graph.plgraph[(*tmp[prid]).vid];
          auto *Arr = vtxA.PE;
          if (vtxA.deg2)
            _mm_prefetch(&graph.pledge[vtxA.offset], _MM_HINT_T0);
          for (uint32 e = 0; e < vtxA.deg1; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              ctx2[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                              (*tmp[prid]).dist);
            else
              ctx2[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*tmp[prid]).dist);
          }
          Arr = graph.pledge;
          for (uint32 e = vtxA.offset; e < vtxA.offset + vtxA.deg2; e += 2) {
            if (!(Arr[e] >> MAX_OFS))
              ctx2[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1],
                                              (*tmp[prid]).dist);
            else
              ctx2[Arr[e] & DW_HALF].emplace_back(
                  MAX_BIT | ((Arr[e] & MIN_HFUP) >> HF_OFST), Arr[e + 1],
                  (*tmp[prid]).dist);
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coro_gather(auto &ctx, auto &low, auto &high) {
    co_yield false;
    for (;;) {
      for (int id = 0; id < (int)low.size(); id += 64) {
        for (int prid = id; prid < std::min(id + 64, (int)low.size()); prid++) {
          _mm_prefetch(&func.vdata[(*low[prid]).e], _MM_HINT_T0);
        }
        co_yield false;
        for (int prid = id; prid < std::min(id + 64, (int)low.size()); prid++) {
          uint32 dst = (*low[prid]).e;
          auto update = func.applyWeight((*low[prid]).w, (*low[prid]).val);
          if (func.gatherFunc(update, dst)) {
            ctx.push_back(func.pushFunc(dst, update));
          }
        }
      }
      for (int id = 0; id < (int)high.size(); id += 64) {
        for (int prid = id; prid < std::min(id + 64, (int)high.size());
             prid++) {
          _mm_prefetch(&graph.highedge[(*high[prid]).w], _MM_HINT_T0);
        }
        co_yield false;
        for (int prid = id; prid < std::min(id + 64, (int)high.size());
             prid++) {
          uint32 deg = (*high[prid]).e & DW_HALF;
          auto dis = (*high[prid]).val;
          for (uint32 e = (*high[prid]).w; e < (*high[prid]).w + deg; e += 2) {
            uint32 dst = graph.highedge[e]; // here
            auto update = func.applyWeight(graph.highedge[e + 1], dis);
            if (func.gatherFunc(update, dst)) {
              ctx.push_back(func.pushFunc(dst, update));
            }
          }
        }
      }
      co_yield true;
    }
  }
  Corobj<bool> coro_sync(auto &ctx2) {
    co_yield false;
    for (;;) {
      wl.scatter(ctx2);
      co_yield true;
    }
  }

public:
  Executor(Func _func, Indexer indexer, Graph &_graph)
      : term(substrate::getSystemTermination(activeThreads)),
        wl(_graph.pnum, indexer), func(_func), graph(_graph) {}

  template <typename Frontier> void initFrontier(Frontier &front) {
    for (auto &a : front) {
      wl.push(a);
    }
    front.clear();
  }

  void print_profiling() {
    wl.summary_cnt();
    wl.summaryTime();
  }

  void initThread() { term.initializeThread(); }

  void operator()() {
    bool isLeader = substrate::ThreadPool::isLeader();
    if (isLeader)
      go<true>();
    else
      go<false>();
  }
};
template <typename OBIM, typename Graph, typename Frontier, typename Indexer,
          typename Func> // RangeTy = LocalRange<InsertBag>
void syncEdgeMap(Graph &graph, Frontier &frontier, const Indexer &indexer,
                 Func func) {
  auto &barrier = getBarrier(activeThreads);
  Executor<OBIM, Func, Indexer, Graph> W(func, indexer, graph);
  W.initFrontier(frontier);
  substrate::getThreadPool().run(activeThreads, [&W]() { W.initThread(); });
  substrate::getThreadPool().run(activeThreads, std::ref(W));
  W.print_profiling();
}
} // end namespace runtime
} // end namespace galois
#endif
