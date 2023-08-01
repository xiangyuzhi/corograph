//
// Created by 15743 on 2023/6/9.
//

#ifndef GALOIS_CONCISE_3_PARTITION_GRAPH_H
#define GALOIS_CONCISE_3_PARTITION_GRAPH_H

#include <fstream>
#include <type_traits>
#include "atomic"
#include "thread"
#include "coroutine"
#include "queue"
#include "unordered_map"
#include "random"
#include "utils.h"

namespace galois::graphs {
    typedef uint32_t uint32;
    typedef uint64_t uint64;
    typedef uint32_t uintPE;
    typedef uint16_t uint16;

#define weighted

#define MAX_NUM 0xffffffff
#define UP_HALF 0xffff0000
#define DW_HALF 0x0000ffff
#define HF_OFST 16
#define MAX_BIT 0x80000000
#define MIN_BIT 0x7fffffff
#define SEC_BIT 0x40000000
#define THD_BIT 0x20000000

#define MIN_BIT 0x7fffffff
#define MIN_HFUP 0x7fff0000
#define MAX_OFS 31
#define SEC_OFS 30
#define THD_OFS 29
#define DW_18bit 0x0003ffff

    inline static void asmPause() {
        asm volatile("pause");
    }

    class SimpleLock {
        mutable std::atomic<int> _lock;
        void slow_lock() const {
            int oldval = 0;
            do {
                while (_lock.load(std::memory_order_acquire) != 0) {
                    asmPause();
                }
                oldval = 0;
            } while (!_lock.compare_exchange_weak(oldval, 1, std::memory_order_acq_rel, std::memory_order_relaxed));
        }

    public:
        constexpr SimpleLock() : _lock(0) {}

        inline void lock() const {
            int oldval = 0;
            if (_lock.load(std::memory_order_relaxed))
                goto slow_path;
            if (!_lock.compare_exchange_weak(oldval, 1, std::memory_order_acq_rel, std::memory_order_relaxed))
                goto slow_path;
            return;
            slow_path:
            slow_lock();
        }

        inline void unlock() const {
            _lock.store(0, std::memory_order_release);
        }

        inline bool try_lock() const {
            int oldval = 0;
            if (_lock.load(std::memory_order_relaxed))
                return false;
            if (!_lock.compare_exchange_weak(oldval, 1, std::memory_order_acq_rel))
                return false;
            return true;
        }

        inline bool is_locked() const {
            return _lock.load(std::memory_order_acquire) & 1;
        }
    };


    typedef std::pair<uint32, uint32> WV;

    class CoroLock{
    public:
        CoroLock() : flag(false) {}
        void lock(std::coroutine_handle<> hd) {
            while (flag.test_and_set(std::memory_order_acquire)) {
                if (!hd.done())
                    hd();
                else std::this_thread::yield();
            }
        }
        void lock() {
            while (flag.test_and_set(std::memory_order_acquire)) {
//            std::this_thread::yield();
                asmPause();
            }
        }
        void unlock() {
            flag.clear(std::memory_order_release);
        }

    private:
        std::atomic_flag flag;
    };


//    template<typename T>
//    class CoroPriorityQueue {
//
//    public:
//        void push(T a, std::coroutine_handle<> hd) {
//            mtx.lock(hd);
//            Q.push(a);
//            mtx.unlock();
//        }
//        void push(T a) {
//            mtx.lock();
//            Q.push(a);
//            mtx.unlock();
//        }
//
//        T pop_top(std::coroutine_handle<> hd){
//            mtx.lock(hd);
//            T ret = Q.top(); //TODO: optimize for new create and only return pointer
//            Q.pop();
//            mtx.unlock();
//            return ret;
//        }
//        T pop_top(){
//            mtx.lock();
//            T ret = Q.top(); //TODO: optimize for new create and only return pointer
//            Q.pop();
//            mtx.unlock();
//            return ret;
//        }
//
//        int size(){
//            return Q.size();
//        }
//
//        T top(){
//            mtx.lock();
//            T ret;
//            if(!Q.size())
//                ret = T(MAX_NUM, MAX_NUM, MAX_NUM, MAX_NUM);
//            else ret = Q.top();
//            mtx.unlock();
//            // here for initial scatter and for empty set to MAX_NUM
//            return ret;
//        }
//
//    private:
//        CoroLock mtx;
//        std::priority_queue<T, std::vector<T>, std::greater<T>> Q;
//    };



    template<typename T>
    class CoroQueue {

    public:
        void lock(){
            mtx.lock();
        }

        void unlock(){
            mtx.unlock();
        }

        void push(T a) {
            mtx.lock();
            Q.push_back(a);
            mtx.unlock();
        }

        void push_lockfree(T a){
            Q.push_back(a);
        }

        T pop_top_lockfree(){
            T ret = Q.front();
            Q.pop_front(); //TODO: optimize for new create and only return pointer
            return ret;
        }

        int size(){
            return Q.size();
        }

    private:
        SimpleLock mtx;
        std::deque<T> Q;
    };

    template<typename T>
    struct part_wrap{//  <e,w,val> | <(1,n),off,val>
        uint32 e,w;
        T val;
    };

    template<typename T>
    struct vertex_warp {
        uint32 vid;
        T dist;
        vertex_warp() {}
        vertex_warp(uint32 _vid, T _dist) : vid(_vid), dist(_dist) {}
        friend bool operator > (const vertex_warp &a, const vertex_warp &b) {
            return a.dist == b.dist ? a.vid > b.vid : a.dist > b.dist;
        }
    };

    template<typename T>
    struct graphPartition{
        uint32 pid; // partition id
//        uint32 numV; // vertex number
        uint32 numE; // edge number

//        uint32* vid; // vertex id list
//        uint32* offset; // CSR offset
        uintPE* ngh; // CSR edge
//        T* weight; // weight

        //CoroQueue<vertex_warp> crq_scatter;// scatter and gather coroqueue
        //CoroQueue<vertex_warp> crq_gather;
        ~graphPartition(){
            free(ngh);
//            free(weight);free(vid);free(offset);
            //free(lh_vid); free(lh_edge);
        }
    };



    struct vtxArr{
        uint16_t deg1, deg2;
        uint32_t PE[14];
        uint32_t offset;
    };// 64B

//    struct vtxArr{
//        uint32 deg1, deg2;
//        uint32 PE[13];
//        uint32 offset;
//    }; // 64B

    template<typename T>
    struct graph{

        uint32 numV; // vertex number
        uint32 numE; // edge number

        uint32* offset; // CSR offset
        uint32* deg; // deg array
        uint32* ngh; // CSR edge
        T* weight; // CSR edge weight

        graphPartition<T>* partition;

//        v2pGraph v2p;
        uint32 psize; // max partition vertex num (1<<16 1<<17 1<<18) [0,1<<16-1] [1<<16, 1<<17-1]
        uint32 plog; // log(partition size)
        uint32 pnum; // partition num

        // power law storage
        vtxArr* plgraph; // powerlaw graph storage
        uint32* pledge;  // powerlaw high graph edge
        uint32* ploffset;  // powerlaw high graph edge
        uint32* plngh;  // powerlaw high graph edge
        uint32 plEnum;
        uint32* highedge;
        uint32 highEnum;
        uint32 plD; // in partition low degree [0,plD]
        uint32 inD; // inline degree, the same as vtxArr arr len

        void init(uint32 _numV, uint32 _numE) {
            numV = _numV;
            numE = _numE;
            offset = new uint32 [_numV+1];
            ngh = new uint32 [_numE];
            deg = new uint32 [_numV];
#ifdef weighted
            weight = new T [_numE];
#endif
        }

        graph() {}
        ~graph(){
            free(offset);
            free(ngh);
//            free(deg);
#ifdef weighted
            free(weight);
#endif
        }

    };
    unsigned int getBinaryLength(unsigned int x) {
        unsigned int length = 0;
        while (x != 0) {
            x >>= 1;
            length++;
        }
        return length;
    }


    template<typename T>
    void load_graph(std::string filename, graph<T> &G) {
        long length;
        char *S = readStringFromFile(filename.c_str(), &length);
        words W = stringToWords(S, length);
        if (strcmp(W.Strings[0], "AdjacencyGraph") != 0) {
            std::cout << "Bad input file: missing header: AdjacencyGraph" << std::endl;
            exit(-1);
        }
        uint64 len = W.m - 1;
        auto *In = new uint32 [len];
        { parallel_for (uint32 i = 0; i < len; i++) In[i] = atol(W.Strings[i + 1]); }
        W.del();
        uint32 n = In[0];
        uint32 m = In[1];

        if (len != n + m + 2) {
            std::cout << "Bad input file: length = " << len << " n+m+2 = " << n + m + 2 << std::endl;
            exit(-1);
        }
        uint32 *offsets = In + 2;
        uint32 *edges = In + 2 + n;

        G.init(n, m);
        printf("Graph numV: %d, numE: %d\n", n, m);

        G.offset[n] = m;
        parallel_for (uint32 i = 0; i < n; i++) {
            uint32 o = offsets[i];
            G.offset[i] = offsets[i];
            uint32 l = ((i == n - 1) ? m : offsets[i + 1]) - offsets[i];
            for (uint32 j = o; j < o + l; j++) {
                G.ngh[j] = edges[j];
            }
        }
        for (uint32 i = 0; i < n; i++) {
            G.deg[i] = G.offset[i + 1] - G.offset[i];
        }
        std::default_random_engine e;
//        std::uniform_int_distribution<int> u(0, getBinaryLength(G.numV));
        std::uniform_int_distribution<int> u(0, 1000);
#ifdef weighted
        for (uint32 j = 0; j < G.numE; j++) {
            G.weight[j] = u(e);
        }
#endif

        uint64 ncnt[14], ecnt[14];
        for (int i = 0; i < 14; i++) ncnt[i] = ecnt[i] = 0;
        for (uint32 i = 0; i < n; i++) {
            if (G.deg[i] < 10) ncnt[G.deg[i]]++, ecnt[G.deg[i]] += G.deg[i];
            else if (G.deg[i] <= 40) ncnt[10]++, ecnt[10] += G.deg[i];
            else ncnt[11]++, ecnt[11] += G.deg[i];
        }
        for (uint32 i = 1; i < 12; i++) {
            ncnt[i] += ncnt[i - 1];
            ecnt[i] += ecnt[i - 1];
        }
        printf("true degree overall\n");
        for (uint32 cc = 0; cc < 10; cc++)
            printf("num: %d cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n", cc, ncnt[cc],
                   ncnt[cc] / (float) G.numV * 100.0, ecnt[cc] / (float) G.numE * 100.0);
        printf("num: 10-40 cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n", ncnt[10],
               ncnt[10] / (float) G.numV * 100.0, ecnt[10] / (float) G.numE * 100.0);
        printf("num: >40 cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n", ncnt[11],
               ncnt[11] / (float) G.numV * 100.0, ecnt[11] / (float) G.numE * 100.0);

        return;
    }


    template<typename T>
    void init_graph(graph<T> &G, commandLine &P) {
        auto filename = P.getOptionValue("-f", "none");
        std::cout << "Reading from file: " << filename << "\n";
        load_graph(filename, G);
        printf("================ Graph load over ================\n");
    }

    
    template<typename T>
    void partition(graph<T> &G, uint32 numthread) {

        uint32 blockVnum = G.numV / (numthread);
        uint32 pow2=1, nlog = 0;
        while(pow2<=blockVnum) {
            nlog ++;
            pow2 <<= 1; 
        }
        nlog --;
        const uint32 plog = 18;
        G.plog = std::min(nlog, plog); printf("plog: %u\n", G.plog);
        uint32 partition_size = 1<<G.plog;
        G.psize = partition_size;
        G.pnum = G.numV / partition_size + 1;
        G.plD = 2; // power law degree
        G.inD = 14; // inplace degree
        printf("psize: %d pnum: %d\n", G.psize, G.pnum);

        G.partition = new graphPartition<T>[G.pnum]();
        graphPartition<T> *P = G.partition;

        uint32* vdeg = new uint32 [G.numV];
        for (uint32 v = 0; v < G.numV; v++) {
            vdeg[v] = 0;
            uint32 prev_e = G.offset[v];
            uint32 prev_p = G.ngh[prev_e] >>G.plog;
            for (uint32 e = G.offset[v]; e < G.offset[v + 1]; e++) {
                uint32 u = G.ngh[e];
                uint32 p = u >>G.plog;
                if (prev_p != p || (e == G.offset[v + 1]-1) ) {
                    if((e == G.offset[v + 1]-1)) e++;
                    uint32 deg = e - prev_e;
                    if(deg <= G.plD){ // low deg
                        vdeg[v] += 2*deg;
                    }
                    else { // high deg
                        vdeg[v] += 2;
                        G.highEnum += 2*deg;
                    }
                    prev_p = p;
                    prev_e = e;
                }
            }
            G.plEnum += vdeg[v];
        }

        for (uint32 i = 0; i < G.pnum; i++) {
            P[i].pid = i;
            P[i].ngh = new uint32 [P[i].numE];
            P[i].numE = 0;
        }
        G.plgraph = (vtxArr*)std::aligned_alloc(64, sizeof(vtxArr)*G.numV);

        uint32 outsum = 0;
        for(uint32 i=0;i<G.numV;i++){
            if(vdeg[i] > G.inD){
                G.plgraph[i].deg1 = G.inD;
                G.plgraph[i].deg2 = vdeg[i] - G.inD;
                G.plgraph[i].offset = outsum;
                outsum += G.plgraph[i].deg2;
            }
            else {
                G.plgraph[i].deg1 = vdeg[i];
                G.plgraph[i].deg2 = 0;
            }
        }
        G.pledge = new uint32[outsum];
//        G.ploffset = new uint32 [G.numV+1];
//        G.plngh = new uint32[G.plEnum];
        G.highedge = new uint32[G.highEnum];

        uint32 outplace_e = 0, high_e = 0;
        for (uint32 v = 0; v < G.numV; v++) {
            uint32 inplace_e = 0 ;
            uint32 prev_e = G.offset[v];
            uint32 prev_p = G.ngh[prev_e] >>G.plog;
            //G.ploffset[v] = in_e;
            for (uint32 e = G.offset[v]; e < G.offset[v + 1]; e++) {
                uint32 u = G.ngh[e];
                uint32 p = u >> G.plog;
                if (prev_p != p  || (e == G.offset[v + 1]-1)) {
                    if((e == G.offset[v + 1]-1)) e++;
                    uint32 deg = e - prev_e;
                    if(deg <= G.plD) { // low deg  // G.plD = 7
                        for (uint32 te = prev_e; te < e; te++) {
                            if(inplace_e<G.inD){ // inline
                                G.plgraph[v].PE[inplace_e ++] = G.ngh[te];
                                G.plgraph[v].PE[inplace_e ++] = G.weight[te];
                            }
                            else { // outline
                                G.pledge[outplace_e++] = G.ngh[te];
                                G.pledge[outplace_e++] = G.weight[te];
                            }
//                            G.plngh[in_e++] = G.ngh[te];
//                            G.plngh[in_e++] = G.weight[te];
                        }
                    }
                    else {
                        if(inplace_e<G.inD){ // inline
                            G.plgraph[v].PE[inplace_e ++] = MAX_BIT | ((deg*2)<<HF_OFST) | prev_p ;
                            G.plgraph[v].PE[inplace_e ++] = high_e;
                        }
                        else { // outline
                            G.pledge[outplace_e++] = MAX_BIT | ((deg*2)<<HF_OFST) | prev_p ;
                            G.pledge[outplace_e++] = high_e;
                        }
                        for (uint32 te = prev_e; te < e; te++) {
                            G.highedge[high_e++] = G.ngh[te];
                            G.highedge[high_e++] = G.weight[te];
                        }
//                        G.plngh[in_e++] = MAX_BIT | ((deg*2)<<HF_OFST) | prev_p ;
//                        G.plngh[in_e++] = out_e;
//                        for (uint32 te = prev_e; te < e; te++) {
//                            G.highedge[out_e++] = G.ngh[te];
//                            G.highedge[out_e++] = G.weight[te];
//                        }
                    }
                    prev_p = p;
                    prev_e = e;
                }
            }
        }
//        G.ploffset[G.numV] = G.plEnum;
//        uint32 out_e = 0;
//        for (uint32 v = 0; v < G.numV; v++) {
//            uint32 in_e = 0;
//            uint32 prev_e = G.offset[v];
//            uint32 prev_p = G.ngh[prev_e] >>G.plog;
//            for (uint32 e = G.offset[v]; e < G.offset[v + 1]; e++) {
//                uint32 u = G.ngh[e];
//                uint32 p = u >> G.plog;
//                if (prev_p != p  || (e == G.offset[v + 1]-1)) {
//                    if((e == G.offset[v + 1]-1)) e++;
//                    uint32 deg = e - prev_e;
//                    if(deg <= G.plD){ // low deg
//                        for(uint32 te = prev_e; te<e; te++){
//                            if(in_e<G.inD){ // inline
//                                G.plgraph[v].PE[in_e ++] = G.ngh[te];
//                                G.plgraph[v].PE[in_e ++] = G.weight[te];
//                            }
//                            else { // outline
//                                G.pledge[out_e++] = G.ngh[te];
//                                G.pledge[out_e++] = G.weight[te];
//                            }
//                        }
//                    }
//                    else { // high partition deg
//                        if(in_e<G.inD){
//                            G.plgraph[v].PE[in_e ++] = prev_p | MAX_BIT;
//                            G.plgraph[v].PE[in_e ++] = P[prev_p].numE;
////                            G.plgraph[v].PE[in_e ++] = G.hPnum;
//                        }
//                        else {
//                            G.pledge[out_e++] = prev_p | MAX_BIT;
//                            G.pledge[out_e++] = P[prev_p].numE;
//                        }
//                        for(uint32 te = prev_e; te<e; te++){
//                            P[prev_p].ngh[P[prev_p].numE++] = G.ngh[te];
//                            P[prev_p].ngh[P[prev_p].numE++] = G.weight[te];
////                            G.hPedge[G.hPnum++] = G.ngh[te];
////                            G.hPedge[G.hPnum++] = G.weight[te];
//                        }
//                        P[prev_p].ngh[P[prev_p].numE-2] |= MAX_BIT;
////                        G.hPedge[G.hPnum-2] |= MAX_BIT;
//                    }
//                    prev_p = p;
//                }
//            }
//        }


        /// profiling distribution
        uint64 ncnt[14],ecnt[14];
        for (int i = 0; i < 14; i++) ncnt[i] = ecnt[i] = 0;
        uint32 sume = 0, bar = 20;
        for (uint32 i = 0; i < G.numV; i++) {
            sume += vdeg[i];
            if ( vdeg[i] < 10) ncnt[vdeg[i]]++, ecnt[vdeg[i]] += vdeg[i];
            else if (vdeg[i] <= bar) ncnt[10]++, ecnt[10] += vdeg[i];
            else ncnt[11]++, ecnt[11] += vdeg[i];
        }
        for (uint32 i = 1; i < 12; i++) {
            ncnt[i] += ncnt[i - 1];
            ecnt[i] += ecnt[i - 1];
        }
        printf("v2p degree overall\n");
        for (uint32 cc = 0; cc < 10; cc++)
            printf("num: %d cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n", cc, ncnt[cc],
                   ncnt[cc] / (float) G.numV * 100.0, ecnt[cc] / (float) sume * 100.0);
        printf("num: 10-%d cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n",bar, ncnt[10],
               ncnt[10] / (float) G.numV * 100.0, ecnt[10] / (float) sume * 100.0);
        printf("num: >%d cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n",bar, ncnt[11],
               ncnt[11] / (float) G.numV * 100.0, ecnt[11] / (float) sume * 100.0);



//        G.partition = new graphPartition<T>[G.pnum]();
////        G.plgraph = new vtxArr[G.numV];
//        G.plgraph = (vtxArr*)std::aligned_alloc(64, sizeof(vtxArr)*G.numV);
//        graphPartition<T> *P = G.partition;
//
//        //uint32* vdeg = new uint32 [G.numV];
//        for (uint32 v = 0; v < G.numV; v++) {
//            //vdeg[v] = 0;
//            uint32 prev_e = G.offset[v];
//            uint32 prev_p = G.ngh[prev_e] >>G.plog;
//            for (uint32 e = G.offset[v]; e < G.offset[v + 1]; e++) {
//                uint32 u = G.ngh[e];
//                uint32 p = u >>G.plog;
//                if (prev_p != p || (e == G.offset[v + 1]-1) ) {
//                    if((e == G.offset[v + 1]-1)) e++;
//                    uint32 deg = e - prev_e;
//                    //vdeg[v] ++;
//                    P[prev_p].numE += 2*deg;
////                    if(deg <= G.plD){ // low deg
////                        vdeg[v] += 2*deg;
////                    }
////                    else { // high deg
////                        vdeg[v] += 2;
////                        P[prev_p].numE += 2*deg;
////                    }
//                    prev_p = p;
//                }
//            }
//        }

//        for (uint32 i = 0; i < G.pnum; i++) {
//            P[i].pid = i;
//            P[i].ngh = new uint32 [P[i].numE];
//            P[i].numE = 0;
//        }

//        uint32 outsum = 0;
//        for(uint32 i=0;i<G.numV;i++){
//            if(G.deg[i] > G.inD){
//                G.plgraph[i].deg1 = G.inD;
//                G.plgraph[i].deg2 = G.deg[i] - G.inD;
//                G.plgraph[i].offset = outsum;
//                outsum += G.plgraph[i].deg2;
//            }
//            else {
//                G.plgraph[i].deg1 = G.deg[i];
//                G.plgraph[i].deg2 = 0;
//            }
////            if(vdeg[i] > G.inD){
////                G.plgraph[i].deg1 = G.inD;
////                G.plgraph[i].deg2 = vdeg[i] - G.inD;
////                G.plgraph[i].offset = outsum;
////                outsum += G.plgraph[i].deg2;
////            }
////            else {
////                G.plgraph[i].deg1 = vdeg[i];
////                G.plgraph[i].deg2 = 0;
////            }
//        }
//        G.pledge = new uint32[outsum];
//        G.pledge = new partition_warp[outsum];

//        uint32 out_e = 0;
//        for (uint32 v = 0; v < G.numV; v++) {
//            uint32 in_e = 0;
//            uint32 prev_e = G.offset[v];
//            uint32 prev_p = G.ngh[prev_e] >>G.plog;
//            for (uint32 e = G.offset[v]; e < G.offset[v + 1]; e++) {
//                uint32 u = G.ngh[e];
//                if (in_e < G.inD) {
//                    G.plgraph[v].PE[in_e++] = u;
//                } else {
//                    G.pledge[out_e++] = u;
//                }
//            }
////                uint32 p = u >> G.plog;
////                if (prev_p != p  || (e == G.offset[v + 1]-1)) {
////                    if((e == G.offset[v + 1]-1)) e++;
////                    if(in_e<G.inD) {
////                        G.plgraph[v].PE[in_e++] = partition_warp(prev_p, prev_e, e);
////                    }
////                    else{
////                        G.pledge[out_e++] = partition_warp(prev_p, prev_e, e);
////                    }
//////                    uint32 st = P[prev_p].numE;
//////                    for(uint32 te = prev_e; te<e; te++){
//////                        P[prev_p].ngh[P[prev_p].numE++] = G.ngh[te];
//////                        P[prev_p].ngh[P[prev_p].numE++] = G.weight[te];
//////                    }
//////                    uint32 ed = P[prev_p].numE;
//////                    if(in_e<G.inD) {
//////                        G.plgraph[v].PE[in_e++] = partition_warp(prev_p, st, ed);
//////                    }
//////                    else{
//////                        G.pledge[out_e++] = partition_warp(prev_p, st, ed);
//////                    }
////                    prev_p = p;
////                }
////            }
//        }


        /// profiling distribution
//        uint64 ncnt[14],ecnt[14];
//        for (int i = 0; i < 14; i++) ncnt[i] = ecnt[i] = 0;
//        uint32 sume = 0;
//        for (uint32 i = 0; i < G.numV; i++) {
//            sume += vdeg[i];
//            if ( vdeg[i] < 10) ncnt[vdeg[i]]++, ecnt[vdeg[i]] += vdeg[i];
//            else if (vdeg[i] <= G.inD) ncnt[10]++, ecnt[10] += vdeg[i];
//            else ncnt[11]++, ecnt[11] += vdeg[i];
//        }
//        for (uint32 i = 1; i < 12; i++) {
//            ncnt[i] += ncnt[i - 1];
//            ecnt[i] += ecnt[i - 1];
//        }
//        printf("v2p degree overall\n");
//        for (uint32 cc = 0; cc < 10; cc++)
//            printf("num: %d cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n", cc, ncnt[cc],
//                   ncnt[cc] / (float) G.numV * 100.0, ecnt[cc] / (float) sume * 100.0);
//        printf("num: 10-%d cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n",G.inD, ncnt[10],
//               ncnt[10] / (float) G.numV * 100.0, ecnt[10] / (float) sume * 100.0);
//        printf("num: >%d cnt: %d vtx rat: %.6f %% edge rat: %.6f %%\n",G.inD, ncnt[11],
//               ncnt[11] / (float) G.numV * 100.0, ecnt[11] / (float) sume * 100.0);

       printf("================ Graph load over ================\n");
    }



}
#endif //GALOIS_CONCISE_3_PARTITION_GRAPH_H
