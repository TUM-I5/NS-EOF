cmake_minimum_required(VERSION 3.21)

set(NSEOF_PROJECT_NAME        "NS-EOF")
set(NSEOF_PROJECT_DESCRIPTION "Navier Stokes - Ernst-Otto-Fischer Teaching Code")
set(NSEOF_AUTHOR_ORGANIZATION "Chair of Scientific Computing, Department of Informatics, Technische Universitaet Muenchen")
set(NSEOF_AUTHOR_DOMAIN       "https://github.com/TUM-I5/NS-EOF")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

project(${NSEOF_PROJECT_NAME}
  DESCRIPTION ${NSEOF_PROJECT_DESCRIPTION}
  HOMEPAGE_URL ${NSEOF_AUTHOR_DOMAIN}
  LANGUAGES C CXX
)

set(CMAKE_C_STANDARD 17 CACHE STRING "")
set(CMAKE_C_STANDARD_REQUIRED TRUE)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "")
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_INLINES_HIDDEN YES)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON CACHE INTERNAL "")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

file(WRITE "${PROJECT_BINARY_DIR}/.gitignore" "*")

include(GNUInstallDirs)
set(DEFAULT_CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

include(TestBigEndian)
test_big_endian(BIG_ENDIAN)

include(ProcessorCount)
ProcessorCount(AVAILABLE_PROCESSOR_COUNT)

include(CTest)
enable_testing()

include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)
if(IPO_SUPPORTED)
  message(STATUS "IPO/LTO enabled")
  set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
else()
  #message(WARNING "IPO/LTO not supported: ${IPO_ERROR}")
  message(WARNING "IPO/LTO not supported")
endif()

add_library(NS-EOF-Interface INTERFACE)
target_compile_features(NS-EOF-Interface INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

target_compile_options(NS-EOF-Interface INTERFACE
  -W
  -Wall
  -Wextra
  -Wpedantic
)

option(ENABLE_SINGLE_PRECISION "Enable single floating-point precision" OFF)
if(ENABLE_SINGLE_PRECISION)
  target_compile_definitions(NS-EOF-Interface INTERFACE ENABLE_SINGLE_PRECISION)
endif()

option(ENABLE_PETSC "Enable the Portable, Extensible Toolkit for Scientific Computation (PETSc)" ON)
if(ENABLE_PETSC)
  find_package(PETSc REQUIRED)
  target_compile_definitions(NS-EOF-Interface INTERFACE ENABLE_PETSC)
  target_link_libraries(NS-EOF-Interface INTERFACE ${PETSc_LIBRARIES})
  target_include_directories(NS-EOF-Interface INTERFACE ${PETSc_INCLUDE_DIRS})
endif()

find_package(MPI REQUIRED)
find_package(Catch2 REQUIRED)
find_package(spdlog REQUIRED)
find_package(tinyxml2 REQUIRED)

target_link_libraries(NS-EOF-Interface INTERFACE
  MPI::MPI_CXX
  spdlog::spdlog
  tinyxml2::tinyxml2
)

add_subdirectory(Source)
add_subdirectory(Tests)
