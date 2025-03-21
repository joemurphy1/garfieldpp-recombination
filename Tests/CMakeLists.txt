add_executable(RandomNumberGenerators RandomNumberGenerators.cpp)
target_link_libraries(RandomNumberGenerators PRIVATE Garfield::Garfield)
add_test(NAME RandomNumberGenerators COMMAND RandomNumberGenerators)

add_executable(Vector.test Vector.test.cpp)
target_link_libraries(Vector.test PRIVATE Garfield::Garfield)
add_test(NAME Vector.test COMMAND Vector.test)

if(GARFIELD_WITH_CUDA)
  add_executable(VectorCuda.test VectorCuda.test.cu)
  target_link_libraries(VectorCuda.test PRIVATE Garfield::Garfield)
  add_test(NAME VectorCuda.test COMMAND VectorCuda.test)
endif()