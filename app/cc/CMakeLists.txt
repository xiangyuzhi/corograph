
add_executable(crg-cc crg-cc.cpp)
target_link_libraries(crg-cc PRIVATE Galois::shmem)


add_executable(crg-cc-async crg-cc-async.cpp)
target_link_libraries(crg-cc-async PRIVATE Galois::shmem)