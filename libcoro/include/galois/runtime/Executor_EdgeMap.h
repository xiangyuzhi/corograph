
#ifndef CRG_EXECUTOR_EDGEMAP_H
#define CRG_EXECUTOR_EDGEMAP_H
#include "galois/Loops.h"
#include "galois/Bag.h"

namespace galois::runtime {

    template<class Graph, class Func>
    class syncExecutor {
        Func func;
        Graph& graph;
        InsertBag<galois::graphs::part_wrap<uint32>> *ctx2;
        InsertBag<uint32> gatherQ;
    public:
        syncExecutor(Graph &_graph, Func _func) : func(_func), graph(_graph) {
            ctx2 = new InsertBag<galois::graphs::part_wrap<uint32>>[graph.pnum];
        }

        template <typename Frontier1, typename Frontier2>
        void EdgeMap(Frontier1& curF, Frontier2& nextF) {
            //printf("scatter...\n");
            struct timespec start, end;
            float time; if( clock_gettime( CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
            galois::do_all(
                galois::iterate(curF),
                [&](const uint32& n) {
                    if(func.filterFunc(n)) return ;
                    auto &vtxA = graph.plgraph[n];
                    auto *Arr = vtxA.PE;
                    if (vtxA.deg2) _mm_prefetch(&graph.pledge[vtxA.offset], _MM_HINT_T0);
                    for (uint32 e = 0; e < vtxA.deg1; e += 2) {
                        if (!(Arr[e] >> MAX_OFS)) ctx2[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1], func.scatterFunc(n));
                        else ctx2[Arr[e] & DW_HALF].emplace_back(MAX_BIT | ((Arr[e] & MIN_HFUP)>>HF_OFST), Arr[e + 1], func.scatterFunc(n));
                    }
                    Arr = graph.pledge;
                    for (uint32 e = vtxA.offset; e < vtxA.offset + vtxA.deg2; e += 2) {
                        if (!(Arr[e] >> MAX_OFS)) ctx2[Arr[e] >> 18].emplace_back(Arr[e], Arr[e + 1], func.scatterFunc(n));
                        else ctx2[Arr[e] & DW_HALF].emplace_back(MAX_BIT | ((Arr[e] & MIN_HFUP)>>HF_OFST), Arr[e + 1], func.scatterFunc(n));
                    }
                },
                galois::steal(), galois::loopname("Merge"));

            for(uint32 i=0;i<graph.pnum;i++){
                if(!ctx2[i].empty()){
                    gatherQ.push(i);
                }
            }
            if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");} time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
            printf("scatter time: %lf sec   ",time);

            if( clock_gettime( CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
            galois::do_all(
                galois::iterate(gatherQ),
                [&](const uint32& n) {
                    for(auto it = ctx2[n].begin(); it!= ctx2[n].end(); it++){
                        if(!((*it).e >> MAX_OFS)){
                            uint32 dst = (*it).e;
                            auto update = func.applyWeight((*it).w, (*it).val);
                            if (func.gatherFunc(update, dst)){
                                nextF.push_back(func.pushFunc(dst, update));
                            }
                        }
                        else {
                            uint32 deg = (*it).e & DW_HALF;
                            auto dis = (*it).val;
                            for (uint32 e = (*it).w; e < (*it).w + deg; e += 2) {
                                uint32 dst = graph.highedge[e];  // here
                                auto update = func.applyWeight(graph.highedge[e + 1], dis);
                                if (func.gatherFunc(update, dst)){
                                    nextF.push_back(func.pushFunc(dst, update));
                                }
                            }
                        }
                    }
                    ctx2[n].clear_serial();
                },
                galois::steal(),galois::chunk_size<1>(),galois::loopname("Merge"));
            gatherQ.clear_serial();
            if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");} time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
            printf("gather time: %lf sec\n",time);
        }

    };




} // namespace galois::runtime

#endif //CRG_EXECUTOR_EDGEMAP_H
