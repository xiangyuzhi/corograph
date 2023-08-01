//
// Created by 15743 on 2023/7/8.
//

#ifndef COROGRAPH_COROBJ_H
#define COROGRAPH_COROBJ_H

#include "coroutine"

template<typename T>
struct Corobj {
    struct promise_type;
    using handle_type = std::coroutine_handle<promise_type>;
    struct promise_type {
        T value_;
        std::exception_ptr exception_;
        Corobj get_return_object() {
            return Corobj(handle_type::from_promise(*this));
        }
        void unhandled_exception() { exception_ = std::current_exception(); }
        std::suspend_never initial_suspend() {return {};}
        std::suspend_always final_suspend() noexcept {return {};}
        template<std::convertible_to<T> From> // C++20 concept
        std::suspend_always yield_value(From &&from) {
            value_ = std::forward<From>(from);
            return {};
        }
        void return_void() {}
    };
    std::coroutine_handle<promise_type> h_;
    Corobj() {}
    Corobj(std::coroutine_handle<promise_type> h) : h_(h) {}
    Corobj& operator=(const Corobj& other) {
        h_ = other.h_;
        return *this;
    }
    explicit operator bool() {
        return !h_.done();
    }
    T operator()() {
        h_();
        return std::move(h_.promise().value_);
    }
};




#endif //COROGRAPH_COROBJ_H
