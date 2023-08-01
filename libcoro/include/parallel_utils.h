//
// Created by 15743 on 2023/2/14.
//

#ifndef EXP_PARALLEL_UTILS_H
#define EXP_PARALLEL_UTILS_H

#include <omp.h>
#define cilk_spawn
#define cilk_sync
#define parallel_main main
#define parallel_for _Pragma("omp parallel for") for
#define parallel_for_1 _Pragma("omp parallel for schedule (static,1)") for
#define parallel_for_256 _Pragma("omp parallel for schedule (static,256)") for

inline void set_num_workers(int n) {
    omp_set_num_threads(n);
}

[[maybe_unused]] static int getWorkers() { return omp_get_max_threads(); }
[[maybe_unused]] static void setWorkers(int n) { omp_set_num_threads(n); }
[[maybe_unused]] static int getWorkerNum() {
    return omp_get_thread_num();
}

#endif //EXP_PARALLEL_UTILS_H
