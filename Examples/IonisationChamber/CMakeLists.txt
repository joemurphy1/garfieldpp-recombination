# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(DriftTube)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

#---Build executable------------------------------------------------------------
add_executable(mdt mdt.C)
target_link_libraries(mdt PRIVATE Garfield::Garfield)

add_executable(Chamber Chamber.C)
target_link_libraries(Chamber PRIVATE Garfield::Garfield)

# FFTW3 library for Poisson solver
find_package(FFTW3 COMPONENTS DOUBLE)
if(FFTW3_FOUND)
  if(DEFINED FFTW3_LIBRARIES)
    target_link_libraries(Chamber PRIVATE ${FFTW3_LIBRARIES})
  else()
    # Try linking by name - often libfftw3 is available via -lfftw3
    target_link_libraries(Chamber PRIVATE fftw3)
  endif()
else()
  # Try the older/alternative FindFFTW
  find_package(FFTW)
  if(FFTW_FOUND)
  target_include_directories(Chamber PRIVATE ${FFTW_INCLUDE_DIRS})
  target_link_libraries(Chamber PRIVATE ${FFTW_LIBRARIES})
  else()
    message(WARNING "FFTW not found. Poisson FFT solver will not be available.")
  endif()
endif()

add_executable(Chamber_vid Chamber_vid.C)
target_link_libraries(Chamber_vid PRIVATE Garfield::Garfield)

add_executable(mdt_mt mdt_mt.C)
target_link_libraries(mdt_mt PRIVATE Garfield::Garfield)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(mdt_mt PRIVATE OpenMP::OpenMP_CXX)
endif()

add_executable(fast_ion_signal fast_ion_signal.C)
target_link_libraries(fast_ion_signal PRIVATE Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file n2_78.08_o2_20.95_ar_0.93_co2_0.04_1atm.gas mdt_elx_delta.txt IonMobility_N2+_N2.txt NegIonMobility_O2-_air.txt)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
