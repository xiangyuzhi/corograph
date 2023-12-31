/*
 * This file belongs to the Galois project, a C++ library for exploiting
 * parallelism. The code is being released under the terms of the 3-Clause BSD
 * License (a copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#ifndef GALOIS_WORKLIST_OBIM_H
#define GALOIS_WORKLIST_OBIM_H

#include <deque>
#include <limits>
#include <type_traits>

#include "galois/FlatMap.h"
#include "galois/runtime/Substrate.h"
#include "galois/substrate/PerThreadStorage.h"
#include "galois/substrate/Termination.h"
#include "galois/worklists/Chunk.h"
#include "galois/worklists/WorkListHelpers.h"
#include "utils.h"

namespace galois {
    namespace worklists {

        template <typename Indexer, typename Container, typename Container2,
                typename Ck, typename Ck2>
        struct OBIM {

            typedef typename Container::value_type value_type;
            typedef typename Indexer::compare Compare;
            typedef typename Container2::value_type part_type;
            typedef uint32 Index;

            Compare compare;
            Index earliest;

            struct ThreadData {
                galois::flat_map<Index, Container*, std::less<Index>> local;
                Index curIndex;
                Index scanStart;
                Container* current;
                unsigned int lastMasterVersion;
                unsigned int numP;
                Container2* currentP;// current owned partition gather ConQueue pointer
                uint32_t Pid;
                Container2** localPQ; // same as the global, for NUMA-aware

                // for profiling
                std::map<std::string, std::map<uint32_t, uint32_t> > cnt;
                std::map<std::string, uint64_t > cntTime;

                ThreadData(Index initial, uint32_t _numP) :
                curIndex(initial), scanStart(initial), current(0), currentP(0), lastMasterVersion(0), numP(_numP){
                    localPQ = new Container2*[_numP];
                }

                void updatePQ(Container2** PQ){
                    for(uint32_t i=0;i<numP; i++){
                        localPQ[i] = PQ[i];
                    }
                }

            };

            substrate::PerThreadStorage<ThreadData> data;
            substrate::PaddedLock<true> masterLock;
            std::deque<std::pair<Index, Container*>> masterLog;

            substrate::PerSocketStorage<ConExtLinkedQueue<Container2 ,true>> gatherQ; // the list of partition id need to be gathered
            Container2** partitionQueue; // per partition queue
            uint32_t numP; // partition num

            std::atomic<unsigned int> masterVersion;
            Indexer indexer;

            bool updateLocal(ThreadData& p) {
                if (p.lastMasterVersion != masterVersion.load(std::memory_order_relaxed)) {
                    for (;p.lastMasterVersion < masterVersion.load(std::memory_order_relaxed); ++p.lastMasterVersion) {
                        std::pair<Index, Container*> logEntry = masterLog[p.lastMasterVersion];
                        p.local[logEntry.first]         = logEntry.second;
                    }
                    return true;
                }
                return false;
            }

            __attribute__((noinline)) Container* slowUpdateLocalOrCreate(ThreadData& p, Index i) {
                do {
                    updateLocal(p);
                    auto it = p.local.find(i);
                    if (it != p.local.end())
                        return it->second;
                } while (!masterLock.try_lock());
                updateLocal(p);
                auto it = p.local.find(i);
                Container* C2 = (it != p.local.end()) ? it->second : nullptr;
                if (!C2) {
                    C2                  = new Container();
                    p.local[i]          = C2;
                    p.lastMasterVersion = masterVersion.load(std::memory_order_relaxed) + 1;
                    masterLog.push_back(std::make_pair(i, C2));
                    masterVersion.fetch_add(1);
                }
                masterLock.unlock();
                return C2;
            }

            __attribute__((noinline)) Ck* slowPop2(ThreadData& p) {
                bool localLeader = substrate::ThreadPool::isLeader();
                Index msS        = this->earliest;
                updateLocal(p);
                msS = p.scanStart;
                if (localLeader) {
                    for (unsigned i = 0; i < runtime::activeThreads; ++i) {
                        Index o = data.getRemote(i)->scanStart;
                        if (this->compare(o, msS))
                            msS = o;
                    }
                } else {
                    Index o = data.getRemote(substrate::ThreadPool::getLeader())->scanStart;
                    if (this->compare(o, msS))
                        msS = o;
                }
                for (auto ii = p.local.lower_bound(msS), ei = p.local.end(); ii != ei; ++ii) {
                    Ck* item = ii->second->pop2();
                    if (item) {
                        p.current   = ii->second;
                        p.curIndex  = ii->first;
                        p.scanStart = ii->first;
                        return item;
                    }
                }
                return nullptr;
            }

            inline Container* updateLocalOrCreate(ThreadData& p, Index i) {
                auto it = p.local.find(i);
                if (it != p.local.end())
                    return it->second;
                return slowUpdateLocalOrCreate(p, i);
            }

        public:
            OBIM(uint32_t _numP, const Indexer& x = Indexer()) : numP(_numP), masterVersion(0), indexer(x), data(this->earliest, _numP),
                                                 earliest(std::numeric_limits<Index>::min()) {
                partitionQueue = new Container2*[_numP];
                for(uint32_t i=0; i<_numP; i++){
                    partitionQueue[i] = new Container2(i);
                }
                for (unsigned i = 0; i < runtime::activeThreads; ++i) {
                    ThreadData &o = *data.getRemote(i);
                    o.updatePQ(partitionQueue);
                }
            }

            ~OBIM() {
                for (auto ii = masterLog.rbegin(), ei = masterLog.rend(); ii != ei; ++ii) {
                    delete ii->second;
                }
            }

            uint32_t getPid(){
                ThreadData &p = *data.getLocal();
                p.Pid = p.currentP->Pid;
                return p.Pid;
            }

            void updateTime(std::string name, auto &st, auto &ed){
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(ed - st);
                ThreadData &p = *data.getLocal();
                p.cntTime[name] += duration.count();
            }
            void summaryTime(){
                uint32_t thd = runtime::activeThreads;
                std::map<std::string, std::vector<uint64> > sumTime;
                for (unsigned i = 0; i < runtime::activeThreads; ++i) {
                    auto& pcnt = data.getRemote(i)->cntTime;
                    for(auto &name : pcnt)
                        sumTime[name.first].push_back(name.second);
                }
                uint64_t totaltime = 0;
                for(auto &name : sumTime){
                    uint64_t tsum = 0;
                    for(auto t : sumTime[name.first]){
                        tsum += t;
                        printf("%lu ",t/1000000);
                    }
                    totaltime += tsum;
                    printf("\n  %s sumtime: %lu ms\n",name.first.c_str(), tsum/1000000);
                }
                printf("totalTime: %lu avgTime: %lu\n", totaltime/1000000, totaltime/1000000/getActiveThreads());
            }

            void update_cnt(std::string name, uint32_t val, uint32_t num = 1){
                ThreadData &p = *data.getLocal();
                p.cnt[name][val] += num;
            }

            void summary_cnt(){
                std::map<std::string, std::map<uint32_t, uint32_t> > sum;
                for (unsigned i = 0; i < runtime::activeThreads; ++i) {
                    auto& pcnt = data.getRemote(i)->cnt;
                    for(auto &name : pcnt)
                        for(auto &val : name.second)
                            sum[name.first][val.first] += val.second;
                }
//                uint32_t level[] = {1000,65536>>1, 65536, 262144>>1, 262144, 262144<<1, 262144<<2, 0xffffffff};
                uint32_t level[] = {10, 100, 1000, 10000, 100000, 1000000};
                uint32_t levelnum = 6;
                for(auto &name : sum){
                    uint32_t cnt[levelnum],scnt[levelnum]; for(uint32_t i=0;i<levelnum;i++) cnt[i]=scnt[i] = 0;
                    printf("summary of %s\n", name.first.c_str());
                    for(auto &it : name.second){
                        for(uint32_t l =0;l<levelnum;l++){
                            if(it.first< level[l]){
                                cnt[l] += it.second;
                                scnt[l] += it.second * it.first;
                                break;
                            }
                        }
                    }
                    uint32_t sumc =0, sums = 0;
                    for(uint32_t i=0;i<levelnum;i++){
                        sumc += cnt[i];
                        sums += scnt[i];
                    }
                    printf("sumc: %llu sums:%llu\n", sumc, sums);
                    for(uint32_t l = 0;l<levelnum;l++){
                        printf("%u-%u: cnt %llu sum %llu %.2f%% %.2f%%\n",
                               (l==0)? 0 : level[l-1], level[l], cnt[l], scnt[l], (double)100*cnt[l]/sumc, (double)100*scnt[l]/sums);
                    }
                }
            }

            void scatter(uint32_t pid, const part_type& pt){
                ThreadData &p = *data.getLocal();
                if(p.localPQ[pid]->push(pt)){
                    gatherQ.getLocal()->push(p.localPQ[pid]);
                }
            }
            void scatter(ThreadData &p, uint32_t pid, const part_type& pt){
                  if(p.localPQ[pid]->push(pt)){ // empty and not in list
                      gatherQ.getLocal()->push(p.localPQ[pid]);
                  }
            }

            template <typename Iter>
            void scatter(ThreadData &p,uint32_t pid, Iter b, Iter e){
                while (b != e)
                    scatter(p, pid, *b++);
            }

            void scatter(auto& ctx2){
                ThreadData &p = *data.getLocal();
                for(uint32_t pid=0;pid<ctx2.size();pid++){
                    //update_cnt("ctx2size", ctx2[i].size());
                    // scatter(p,i,ctx2[i].begin(), ctx2[i].end());
                    for(auto iter = ctx2[pid].begin(); iter != ctx2[pid].end();){
                        if(p.localPQ[pid]->push(*iter++)){ // empty and not in list
                            gatherQ.getLocal()->push(p.localPQ[pid]);
                        }
                    }
                    ctx2[pid].clear();
                }
            }

            void sync(){
                ThreadData &p = *data.getLocal();
                for(uint32_t i=0;i<numP;i++){
                    if(p.localPQ[i]->pushNext()){ // push all local update to partition queue
                        gatherQ.getLocal()->push(p.localPQ[i]);
                    }
                }
            }

            Container2* popPartByID(unsigned int i) {
                auto& I = *gatherQ.getRemote(i);
                return I.pop();
            }

            Container2* popPart() {
                int id   = substrate::ThreadPool::getTID(); // thread id
                Container2* r = popPartByID(id);
                if (r)
                    return r;

                for (int i = id + 1; i < (int)gatherQ.size(); ++i) {
                    r = popPartByID(i);
                    if (r)
                        return r;
                }
                for (int i = 0; i < id; ++i) {
                    r = popPartByID(i);
                    if (r)
                        return r;
                }
                return 0;
            }


            bool pop_part(){ // try to pop a partition from
                ThreadData& p = *data.getLocal();
                p.currentP = popPart();//gatherQ.getLocal()->pop();
                if( p.currentP ){
                    return true;
                }
                return false;
            }

            Ck2* pop_gather(){
                ThreadData& p = *data.getLocal();
                Ck2* item = p.currentP->pop();
                return item;
            }

            void push(const value_type& val) {
                Index index   = indexer(val);
                ThreadData& p = *data.getLocal();

                if (index == p.curIndex && p.current) {
                    p.current->push(val);
                    return;
                }

                Container* C = updateLocalOrCreate(p, index);
                if (this->compare(index, p.scanStart))
                    p.scanStart = index;
                if (this->compare(index, p.curIndex)) {
                    p.curIndex = index;
                    p.current  = C;
                }
                C->push(val);
            }

            template <typename Iter>
            void push(Iter b, Iter e) {
                while (b != e)
                    push(*b++);
            }

            Ck* pop2(){ // pop a chunk pointer == pop3 + slowPop3
                ThreadData& p = *data.getLocal();
                Container* C        = p.current;
                Ck* item = nullptr;
                if(C) item = C->pop2();
                if(item) return item;
                return slowPop2(p);
            }

            Ck* pop3(){ // pop a chunk pointer
                ThreadData& p = *data.getLocal();
                Container* C        = p.current;
                Ck* item = nullptr;
                if(C) item = C->pop2();
                return item;
            }

            __attribute__((noinline)) Ck* slowPop3() {
                ThreadData& p = *data.getLocal();
                bool localLeader = substrate::ThreadPool::isLeader();
                Index msS        = this->earliest;
                updateLocal(p);
                msS = p.scanStart;
                if (localLeader) {
                    for (unsigned i = 0; i < runtime::activeThreads; ++i) {
                        Index o = data.getRemote(i)->scanStart;
                        if (this->compare(o, msS))
                            msS = o;
                    }
                } else {
                    Index o = data.getRemote(substrate::ThreadPool::getLeader())->scanStart;
                    if (this->compare(o, msS))
                        msS = o;
                }
                for (auto ii = p.local.lower_bound(msS), ei = p.local.end(); ii != ei; ++ii) {
                    Ck* item = ii->second->pop2();
                    if (item) {
                        p.current   = ii->second;
                        p.curIndex  = ii->first;
                        p.scanStart = ii->first;
                        return item;
                    }
                }
                return nullptr;
            }
        };
    } // end namespace worklists
} // end namespace galois

#endif
