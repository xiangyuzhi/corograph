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

#ifndef GALOIS_WORKLIST_CHUNK_H
#define GALOIS_WORKLIST_CHUNK_H

#include "galois/config.h"
#include "galois/FixedSizeRing.h"
#include "galois/runtime/Mem.h"
#include "galois/substrate/PaddedLock.h"
#include "galois/worklists/WLCompileCheck.h"
#include "galois/worklists/WorkListHelpers.h"

namespace galois {
    namespace runtime {
        extern unsigned activeThreads;
    }
    namespace worklists {

        namespace internal {

            template <typename T, int ChunkSize>
            class MyChunk : public FixedSizeRing2<T, ChunkSize> {
                MyChunk* next;
            public:
                MyChunk() : next(0) {}
                MyChunk*& getNext() { return next; }
                MyChunk* const& getNext() const { return next; }
            };


            template <typename T, int ChunkSize>
            struct ChunkMaster2 :  public ConExtLinkedQueue<ChunkMaster2<T,ChunkSize> ,true>::ListNode{
            public:
                using Chunk = MyChunk<T, ChunkSize>;
                runtime::FixedSizeAllocator<Chunk> alloc;
                struct p {
                    Chunk* cur;
                    Chunk* next;
                    p() : cur(0), next(0) {}
                };
                substrate::PerThreadStorage<p> data;
                ConExtLinkedQueue2<Chunk, true> Q;
                uint32_t Pid;

                Chunk* mkChunk() {
                    Chunk* ptr = alloc.allocate(1);
                    alloc.construct(ptr);
                    return ptr;
                }

                void delChunk(Chunk* ptr) {
                    alloc.destroy(ptr);
                    alloc.deallocate(ptr, 1);
                }
                bool pushChunk(Chunk* C) {
                    return Q.push(C);
                }
                Chunk* popChunk() {
                    return Q.pop();
                }

                template <typename... Args>
                bool emplacei(p& n, Args&&... args) {
                    if (n.next && !n.next->emplace_back(std::forward<Args>(args)...))
                        return false;
                    if (n.next){ // full
                        bool ret = pushChunk(n.next);
                        n.next = 0;
                        return ret;
                    }
                    n.next = mkChunk();
                    n.next->emplace_back(std::forward<Args>(args)...);
                    return false;
                }

                typedef T value_type;

                explicit ChunkMaster2(uint32_t _pid) : Pid(_pid) {}
                ChunkMaster2(const ChunkMaster2&) = delete;
                ChunkMaster2& operator=(const ChunkMaster2&) = delete;

                bool push(const value_type& val) {
                    p& n = *data.getLocal();
                    return emplacei(n, val);
                }

                bool pushNext(){
                    bool ret = false;
                    p& n = *data.getLocal();
                    if (n.next){//&& !n.next->empty()
                        ret = pushChunk(n.next);
                        n.next = 0;
                        //n.next = mkChunk();
                    }
                    return ret;
                }

                Chunk* pop() {
                    p& n = *data.getLocal();
//                    if (n.cur && !n.cur->empty()){
//                        return n.cur;
//                    }
                    if (n.cur)
                        delChunk(n.cur);
                    n.cur = popChunk();
                    if (!n.cur) {
                        n.cur  = n.next;
                        n.next = 0;
                    }
                    return n.cur;
//                    if (n.cur && !n.cur->empty())
//                        return n.cur;
//                    return nullptr;
                }
            };

            template <typename T, int ChunkSize>
            struct ChunkMaster {
            public:
                using Chunk = MyChunk<T,ChunkSize>;
                runtime::FixedSizeAllocator<Chunk> alloc;

                struct p {
                    Chunk* cur;
                    Chunk* next;
                    p() : cur(0), next(0) {}
                };

                substrate::PerThreadStorage<p> data;
                substrate::PerSocketStorage<ConExtLinkedQueue<Chunk, true>> Q;

                Chunk* mkChunk() {
                    Chunk* ptr = alloc.allocate(1);
                    alloc.construct(ptr);
                    return ptr;
                }

                void delChunk(Chunk* ptr) {
                    alloc.destroy(ptr);
                    alloc.deallocate(ptr, 1);
                }

                void pushChunk(Chunk* C) {
                    auto& I = *Q.getLocal();
                    I.push(C);
                }

                Chunk* popChunkByID(unsigned int i) {
                    auto& I = *Q.getRemote(i);
                    return I.pop();
                }

                Chunk* popChunk() {
                    int id   = substrate::ThreadPool::getTID(); // thread id
                    Chunk* r = popChunkByID(id);
                    if (r)
                        return r;

                    for (int i = id + 1; i < (int)Q.size(); ++i) {
                        r = popChunkByID(i);
                        if (r)
                            return r;
                    }
                    for (int i = 0; i < id; ++i) {
                        r = popChunkByID(i);
                        if (r)
                            return r;
                    }
                    return 0;
                }

                template <typename... Args>
                void emplacei(p& n, Args&&... args) {
                    if (n.next && !n.next->emplace_back(std::forward<Args>(args)...))
                        return ;
                    if (n.next){
                        pushChunk(n.next);
                        n.next = 0;
                        return ;
                    }
                    n.next = mkChunk();
                    n.next->emplace_back(std::forward<Args>(args)...);
                }


                typedef T value_type;
                typedef MyChunk<T, ChunkSize> Ck;

                ChunkMaster()                   = default;
                ChunkMaster(const ChunkMaster&) = delete;
                ChunkMaster& operator=(const ChunkMaster&) = delete;

                void push(const value_type& val) {
                    p& n = *data.getLocal();
                    emplacei(n, val);
                }

                template <typename Iter>
                void push(Iter b, Iter e) {
                    p& n = *data.getLocal();
                    while (b != e)
                        emplacei(n, *b++);
                }

                Ck* pop2() {
                    p& n = *data.getLocal();
//                    if (n.cur && !n.cur->empty()){
//                        return n.cur;
//                    }
                    if (n.cur)
                        delChunk(n.cur);
                    n.cur = popChunk();
                    if (!n.cur) {
                        n.cur  = n.next;
                        n.next = 0;
                    }
                    return n.cur;
//                    if (n.cur && !n.cur->empty())
//                        return n.cur;
//                    return nullptr;
                }

//                galois::optional<value_type> pop() {
//                    p& n = *data.getLocal();
//                    galois::optional<value_type> retval;
//                    if (n.cur && (retval = n.cur->extract_front()))
//                        return retval;
//                    if (n.cur)
//                        delChunk(n.cur);
//                    n.cur = popChunk();
//                    if (!n.cur) {
//                        n.cur  = n.next;
//                        n.next = 0;
//                    }
//                    if (n.cur)
//                        return n.cur->extract_front();
//                    return galois::optional<value_type>();
//                }
            };

        } // namespace internal
        template <int ChunkSize = 64, typename T = int> using CM = internal::ChunkMaster<T, ChunkSize>;
        template <int ChunkSize = 64, typename T = int> using CM2 = internal::ChunkMaster2<T, ChunkSize>;
        template <int ChunkSize = 64, typename T = int> using CK = internal::MyChunk<T, ChunkSize>;
    } // end namespace worklists
} // end namespace galois

#endif
