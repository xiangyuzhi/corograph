
add_executable(crg-kcore-single crg-kcore-single.cpp)
target_link_libraries(crg-kcore-single PRIVATE Galois::shmem)

 add_executable(crg-kcore crg-kcore.cpp)
 target_link_libraries(crg-kcore PRIVATE Galois::shmem)