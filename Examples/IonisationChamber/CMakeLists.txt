# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(IonisationChamber)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

#---Build executable------------------------------------------------------------
add_executable(plots plots.C)
target_link_libraries(plots PRIVATE Garfield::Garfield)

add_executable(Chamber Chamber.C)
target_link_libraries(Chamber PRIVATE Garfield::Garfield)

add_executable(RoosChamber RoosChamber.C)
target_link_libraries(RoosChamber PRIVATE Garfield::Garfield)

# FFTW3 library for Poisson solver
find_package(FFTW3 COMPONENTS DOUBLE)
if(FFTW3_FOUND)
  if(DEFINED FFTW3_LIBRARIES)
    target_link_libraries(Chamber PRIVATE ${FFTW3_LIBRARIES})
    target_link_libraries(RoosChamber PRIVATE ${FFTW3_LIBRARIES})
  else()
    target_link_libraries(Chamber PRIVATE fftw3)
    target_link_libraries(RoosChamber PRIVATE fftw3)
  endif()
  find_package(Threads REQUIRED)
  target_link_libraries(Chamber PRIVATE fftw3_threads Threads::Threads)
  target_link_libraries(RoosChamber PRIVATE fftw3_threads Threads::Threads)
else()
  find_package(FFTW)
  if(FFTW_FOUND)
    target_include_directories(Chamber PRIVATE ${FFTW_INCLUDE_DIRS})
    target_link_libraries(Chamber PRIVATE ${FFTW_LIBRARIES})
    target_include_directories(RoosChamber PRIVATE ${FFTW_INCLUDE_DIRS})
    target_link_libraries(RoosChamber PRIVATE ${FFTW_LIBRARIES})
  else()
    message(WARNING "FFTW not found. Poisson FFT solver will not be available.")
  endif()
endif()

add_executable(Chamber_vid Chamber_vid.C)
target_link_libraries(Chamber_vid PRIVATE Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file IonMobility_N2+_N2.txt NegIonMobility_O2-_air.txt)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

# Copy all gas files from the gas folder
file(GLOB _gas_files "${CMAKE_CURRENT_SOURCE_DIR}/gas_files/*.gas")
foreach(_file ${_gas_files})
  get_filename_component(_filename ${_file} NAME)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/gas_files/${_filename} COPYONLY)
endforeach()

# Copy all files from the Elmer folder
file(GLOB _elmer_files "${CMAKE_CURRENT_SOURCE_DIR}/Elmer/*")
foreach(_file ${_elmer_files})
  get_filename_component(_filename ${_file} NAME)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/Elmer/${_filename} COPYONLY)
endforeach()