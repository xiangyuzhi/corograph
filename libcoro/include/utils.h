//
// Created by 15743 on 2023/2/14.
//
#ifndef UTILS_H
#define UTILS_H
#include "parallel_utils.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include "vector"

typedef uint32_t uintPE;
typedef uint32_t uint32;
typedef uint16_t uint16;
typedef uint64_t uint64;

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

#define newA(__E, __n) (__E*) malloc((__n)*sizeof(__E))


struct words {
    char *Chars;  // array storing all strings
    long n; // total number of characters
    char **Strings; // pointers to strings (all should be null terminated)
    long m; // number of substrings
    words() {}

    words(char *C, long nn, char **S, long mm)
            : Chars(C), n(nn), Strings(S), m(mm) {}

    void del() {
        free(Chars);
        free(Strings);
    }
};



struct commandLine {
    int argc;
    char **argv;
    std::string comLine;

    commandLine(int _c, char **_v, std::string _cl)
            : argc(_c), argv(_v), comLine(_cl) {
        if (getOption("-h") || getOption("-help"))
            badArgument();
    }

    commandLine(int _c, char **_v)
            : argc(_c), argv(_v), comLine("bad arguments") {}

    void badArgument() {
        std::cout << "usage: " << argv[0] << " " << comLine << std::endl;
        exit(0);
    }

    // get an argument
    // i is indexed from the last argument = 0, second to last indexed 1, ..
    char *getArgument(int i) {
        if (argc < 2 + i) badArgument();
        return argv[argc - 1 - i];
    }

    // looks for two filenames
    std::pair<char *, char *> IOFileNames() {
        if (argc < 3) badArgument();
        return std::pair<char *, char *>(argv[argc - 2], argv[argc - 1]);
    }

    std::pair<size_t, char *> sizeAndFileName() {
        if (argc < 3) badArgument();
        return std::pair<size_t, char *>(std::atoi(argv[argc - 2]), (char *) argv[argc - 1]);
    }

    bool getOption(std::string option) {
        for (int i = 1; i < argc; i++)
            if ((std::string) argv[i] == option) return true;
        return false;
    }

    char *getOptionValue(std::string option) {
        for (int i = 1; i < argc - 1; i++)
            if ((std::string) argv[i] == option) return argv[i + 1];
        return NULL;
    }

    std::string getOptionValue(std::string option, std::string defaultValue) {
        for (int i = 1; i < argc - 1; i++)
            if ((std::string) argv[i] == option) return (std::string) argv[i + 1];
        return defaultValue;
    }

    long getOptionLongValue(std::string option, long defaultValue) {
        for (int i = 1; i < argc - 1; i++)
            if ((std::string) argv[i] == option) {
                long r = atol(argv[i + 1]);
                if (r < 0) badArgument();
                return r;
            }
        return defaultValue;
    }

    int getOptionIntValue(std::string option, int defaultValue) {
        for (int i = 1; i < argc - 1; i++)
            if ((std::string) argv[i] == option) {
                int r = atoi(argv[i + 1]);
                if (r < 0) badArgument();
                return r;
            }
        return defaultValue;
    }

    double getOptionDoubleValue(std::string option, double defaultValue) {
        for (int i = 1; i < argc - 1; i++)
            if ((std::string) argv[i] == option) {
                double val;
                if (sscanf(argv[i + 1], "%lf", &val) == EOF) {
                    badArgument();
                }
                return val;
            }
        return defaultValue;
    }

};


template <class ET>
inline bool CAS(ET *ptr, ET oldv, ET newv) {
    if (sizeof(ET) == 1) {
        return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv), *((bool*)&newv));
    } else if (sizeof(ET) == 4) {
        return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv), *((int*)&newv));
    } else if (sizeof(ET) == 8) {
        return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv), *((long*)&newv));
    }
    else {
        std::cout << "CAS bad length : " << sizeof(ET) << std::endl;
        abort();
    }
}

template <class ET>
inline bool writeMin(ET *a, ET b) {
    ET c; bool r=0;
    do c = *a;
    while (c > b && !(r=CAS(a,c,b)));
    return r;
}


inline bool isSpace(char c);

unsigned int getBinaryLength(unsigned int x);
words stringToWords(char *Str, uint64_t n);

char *readStringFromFile(const char *fileName, long *length);


#endif //UTILS_H
