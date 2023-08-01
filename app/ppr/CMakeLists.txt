add_executable(crg-ppr crg-ppr.cpp)
target_link_libraries(crg-ppr PRIVATE Galois::shmem)