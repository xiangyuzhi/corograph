
add_executable(crg-sssp crg-sssp.cpp)
target_link_libraries(crg-sssp PRIVATE Galois::shmem)