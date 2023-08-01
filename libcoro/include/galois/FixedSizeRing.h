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

#ifndef GALOIS_FIXEDSIZERING_H
#define GALOIS_FIXEDSIZERING_H

#include <atomic>
#include <utility>

#include <boost/mpl/if.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/reverse_iterator.hpp>

#include "galois/config.h"
#include "galois/optional.h"
#include "galois/LazyArray.h"

namespace galois {

template <typename T, unsigned ChunkSize = 64>
class FixedSizeRing2 {
public:
    LazyArray<T, ChunkSize> datac;
    unsigned count;
    T* at(unsigned i) { return &datac[i]; }
    const T* at(unsigned i) const { return &datac[i]; }
public:
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
    typedef const T& const_reference;

    FixedSizeRing2() : count(0) {}
    ~FixedSizeRing2() { clear(); }
//    bool empty() const {
//        return count == 0;
//    }
    bool full() const {
        return count == ChunkSize;
    }
    reference getAt(unsigned x) {
        return *at(x);
    }
    const_reference getAt(unsigned x) const {
        return *at(x);
    }
    void clear() {
        for (unsigned x = 0; x < count; ++x)
            datac.destroy(x);
        count = 0;
    }
    template <typename U>
    pointer push_back(U&& val) {
        return emplace_back(std::forward<U>(val));
    }
//    template <typename... Args>
//    pointer emplace_back(Args&&... args) {
//        if (full())
//            return 0;
//        return datac.emplace(count++, std::forward<Args>(args)...);
//    }
    template <typename... Args>
    bool emplace_back(Args&&... args) {
        datac.emplace(count++, std::forward<Args>(args)...);
        return full();
    }

    galois::optional<value_type> extract_at(int x) {
        return galois::optional<value_type>(*at(x));
    }
};


//! Ordered collection of bounded size
template <typename T, unsigned ChunkSize = 64>
class FixedSizeRing {
public:
  LazyArray<T, ChunkSize> datac;
  unsigned start;
  unsigned count;

  T* at(unsigned i) { return &datac[i]; }
  const T* at(unsigned i) const { return &datac[i]; }

  bool precondition() const { return count <= ChunkSize && start <= ChunkSize; }

  template <typename U>
  class Iterator
      : public boost::iterator_facade<Iterator<U>, U,
                                      boost::random_access_traversal_tag> {
    friend class boost::iterator_core_access;
    U* base;
    unsigned cur;
    unsigned count;

    template <typename OtherTy>
    bool equal(const Iterator<OtherTy>& o) const {
//      assert(base && o.base);
      return &base[cur] == &o.base[o.cur] && count == o.count;
    }

    U& dereference() const { return base[cur]; }

    void increment() {
//      assert(base && count != 0);
      count -= 1;
      cur = (cur + 1) % ChunkSize;
    }

    void decrement() {
//      assert(base && count < ChunkSize);
      count += 1;
      cur = (cur + ChunkSize - 1) % ChunkSize;
    }

    void advance(ptrdiff_t x) {
      count -= x;
      cur = (cur + ChunkSize + x) % ChunkSize;
    }

    ptrdiff_t distance_to(const Iterator& o) const {
      ptrdiff_t c  = count;
      ptrdiff_t oc = o.count;
      return c - oc;
    }


    Iterator() : base(0), cur(0), count(0) {}

    template <typename OtherTy>
    Iterator(const Iterator<OtherTy>& o)
        : base(o.base), cur(o.cur), count(o.count) {}

    Iterator(U* b, unsigned c, unsigned co) : base(b), cur(c), count(co) {}
  };

public:
  typedef T value_type;
  typedef T* pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef Iterator<T> iterator;
  typedef Iterator<const T> const_iterator;
  typedef boost::reverse_iterator<Iterator<T>> reverse_iterator;
  typedef boost::reverse_iterator<Iterator<const T>> const_reverse_iterator;

  FixedSizeRing() : start(0), count(0) {}

  template <typename InputIterator>
  FixedSizeRing(InputIterator first, InputIterator last) : start(0), count(0) {
    while (first != last) {
//      assert(count < ChunkSize);
      datac.emplace(count++, *first++);
    }
  }

  FixedSizeRing(const FixedSizeRing& o) = delete;
  FixedSizeRing& operator=(const FixedSizeRing& o) = delete;

  ~FixedSizeRing() { clear(); }

  unsigned size() const {
//    assert(precondition());
    return count;
  }

  bool empty() const {
    //assert(precondition());
    return count == 0;
  }

  bool full() const {
//    assert(precondition());
    return count == ChunkSize;
  }

  reference getAt(unsigned x) {
//    assert(precondition());
//    assert(!empty());
    return *at(start + x);
  }

  const_reference getAt(unsigned x) const {
//    assert(precondition());
//    assert(!empty());
    return *at(start + x);
  }

  void clear() {
//    assert(precondition());
    for (unsigned x = 0; x < count; ++x)
      datac.destroy(start + x);
    count = 0;
    start = 0;
  }

  // NB(ddn): Keeping emplace_front/_back code paths separate to improve
  // branch prediction etc
  template <typename... Args>
  pointer emplace(iterator pos, Args&&... args) {
    if (full())
      return 0;
    unsigned i;
    if (pos == begin()) {
      i = start = (start + ChunkSize - 1) % ChunkSize;
      ++count;
    } else if (pos == end()) {
      i = (start + count) % ChunkSize;
      ++count;
    } else {
      auto d = std::distance(begin(), pos);
      i      = (start + d) % ChunkSize;
      emplace_back();
      std::move_backward(begin() + d, end() - 1, end());
      datac.destroy(i);
    }
    return datac.emplace(i, std::forward<Args>(args)...);
  }

  template <typename U>
  pointer push_front(U&& val) {
    return emplace_front(std::forward<U>(val));
  }

  template <typename... Args>
  pointer emplace_front(Args&&... args) {
    if (full())
      return 0;
    start = (start + ChunkSize - 1) % ChunkSize;
    ++count;
    return datac.emplace(start, std::forward<Args>(args)...);
  }

  template <typename U>
  pointer push_back(U&& val) {
    return emplace_back(std::forward<U>(val));
  }

  template <typename... Args>
  pointer emplace_back(Args&&... args) {
    if (full())
      return 0;
    unsigned end = start + count;//(start + count) % ChunkSize;
    ++count;
    return datac.emplace(end, std::forward<Args>(args)...);
  }

  reference front() {
//    assert(precondition());
//    assert(!empty());
    return *at(start);
  }

  const_reference front() const {
//    assert(precondition());
//    assert(!empty());
    return *at(start);
  }

    galois::optional<value_type> extract_at(int x) {
        return galois::optional<value_type>(*at(x));
    }

//  galois::optional<value_type> extract_front() {
//    if (!empty()) {
//      galois::optional<value_type> retval(front());
//      pop_front();
//      return retval;
//    }
//    return galois::optional<value_type>();
//  }

  void pop_front() {
//    assert(precondition());
//    assert(!empty());
    datac.destroy(start++);
//    start = (start + 1) % ChunkSize;
    --count;
  }

  reference back() {
//    assert(precondition());
//    assert(!empty());
    return *at((start + count - 1) % ChunkSize);
  }

  const_reference back() const {
//    assert(precondition());
//    assert(!empty());
    return *at((start + count - 1) % ChunkSize);
  }

  galois::optional<value_type> extract_back() {
    if (!empty()) {
      galois::optional<value_type> retval(back());
      pop_back();
      return retval;
    }
    return galois::optional<value_type>();
  }

  void pop_back() {
//    assert(precondition());
//    assert(!empty());
    datac.destroy((start + count - 1) % ChunkSize);
    --count;
  }

  iterator begin() { return iterator(at(0), start, count); }
  iterator end() { return iterator(at(0), (start + count) % ChunkSize, 0); }
  const_iterator begin() const { return const_iterator(at(0), start, count); }
  const_iterator end() const {
    return const_iterator(at(0), (start + count) % ChunkSize, 0);
  }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_iterator rbegin() const { const_reverse_iterator(this->end()); }
  const_iterator rend() const { const_reverse_iterator(this->begin()); }
};

} // namespace galois
#endif
