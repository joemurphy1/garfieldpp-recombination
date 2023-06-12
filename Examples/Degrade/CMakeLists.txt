cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Degrade)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

# ---Build executable------------------------------------------------------------
add_executable(mip mip.C)
target_link_libraries(mip Garfield::Garfield)
