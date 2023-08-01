# Overview
CoroGraph: Bridging Work Efficiency and Cache Efficiency for Graph Algorithms


CoroGraph is a graph framework implemented based on [Galois](https://github.com/IntelligentSoftwareSystems/Galois).
CoroGraph aims to optimize the cache locality for irregular access pattern while maintain the graph preferred
execution model (work effciency).

# Building CoroGraph

## Dependencies

CoroGraph builds, runs, and has been tested on GNU/Linux.

- A modern C++ compiler compliant with the C++-20 standard (gcc >= 11)
- CMake (>= 3.13)
- Boost library (>= 1.58.0, we recommend building/installing the full library)
- libfmt (>= 4.0)

## Compiling and Running CoroGraph


```
SRC_DIR=`pwd` # Or top-level Galois source dir
BUILD_DIR=<path-to-your-build-dir>

mkdir -p $BUILD_DIR
cmake -S $SRC_DIR -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Release
cd $BUILD_DIR
make ..
```

## Graph Input Format
We use the adjacent graph format from the Problem Based Benchmark Suite (http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html).

## Running SSSP as Example

```
cd $BUILD_DIR
./app/sssp -f /path/to/your/graph.adj -delta <delta> -t <threads num>
```





