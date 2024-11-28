cmake_minimum_required(VERSION 3.27)
project(OccultationSimulator VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(EXECUTABLE_NAME slc)
set(LIBRARY_NAME occultation_simulator)

add_subdirectory(occultation_simulator)
add_subdirectory(cli)
add_subdirectory(external)

include(CTest)
add_subdirectory(tests)
