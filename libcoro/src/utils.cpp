//
// Created by 15743 on 2023/6/10.
//
#include "utils.h"


uint32_t rand_in_range(uint32_t max) { return rand() % max; }


// A structure that keeps a sequence of strings all allocated from
// the same block of memory

inline bool isSpace(char c) {
    switch (c) {
        case '\r':
        case '\t':
        case '\n':
        case 0:
        case ' ' :
            return true;
        default :
            return false;
    }
}
unsigned int getBinaryLength(unsigned int x) {
    unsigned int length = 0;
    while (x != 0) {
        x >>= 1;
        length++;
    }
    return length;
}

words stringToWords(char *Str, uint64_t n) {
    parallel_for (uint64_t i = 0; i < n; i++)if (isSpace(Str[i])) Str[i] = 0;

    // mark start of words
    bool *FL = newA(bool, n);
    FL[0] = Str[0];
    parallel_for (uint64_t i = 1; i < n; i++) FL[i] = Str[i] && !Str[i - 1];

    uint32_t worker_count = getWorkers();
    std::vector<uint64_t> sub_counts(worker_count, 0);
    uint64_t section_count = (n / worker_count) + 1;
    parallel_for_1 (uint64_t i = 0; i < worker_count; i++) {
        uint64_t start = i * section_count;
        uint64_t end = std::min((i + 1) * section_count, n);
        uint64_t local_count = 0;
        for (uint64_t j = start; j < end; j++) {
            if (FL[j]) {
                local_count += 1;
            }
        }
        sub_counts[i] = local_count;
    }
    // count and prefix sum
    for (uint32_t i = 1; i < worker_count; i++) {
        sub_counts[i] += sub_counts[i - 1];
    }
    uint64_t m = sub_counts[worker_count - 1];
    uint64_t *offsets = newA(uint64_t, m);
    parallel_for_1 (uint64_t i = 0; i < worker_count; i++) {
        uint64_t start = i * section_count;
        uint64_t end = std::min((i + 1) * section_count, n);
        uint64_t offset;
        if (i == 0) offset = 0;
        else offset = sub_counts[i - 1];
        for (uint64_t j = start; j < end; j++) {
            if (FL[j] == 1) {
                offsets[offset++] = j;
            }
        }
    }

    // pointer to each start of word
    char **SA = newA(char*, m);
    parallel_for (uint64_t j = 0; j < m; j++) SA[j] = Str + offsets[j];

    free(offsets);
    free(FL);
    return words(Str, n, SA, m);
}


char *readStringFromFile(const char *fileName, long *length) {
    std::ifstream file(fileName, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        std::cout << "Unable to open file: " << fileName << std::endl;
        abort();
    }
    long end = file.tellg();
    file.seekg(0, std::ios::beg);
    long n = end - file.tellg();
    char *bytes = (char *) malloc((n + 1) * sizeof(char));
    file.read(bytes, n);
    file.close();
    *length = n;
    return bytes;
}


