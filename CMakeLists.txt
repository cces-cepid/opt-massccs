cmake_minimum_required(VERSION 3.10)

project(massccs
    VERSION         2.0
    DESCRIPTION     "massccs ompc version"
    LANGUAGES       C CXX)

set(CMAKE_C_STANDARD   99)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Remove annoying warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-return-type -Wno-absolute-value")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type -Wno-absolute-value")

# Ensure we are using OMPC Clang
# ------------------------------------------------------------------------------

if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
    execute_process(
        COMMAND ${CMAKE_C_COMPILER} --version
        OUTPUT_VARIABLE CLANG_VERSION)
    if(CLANG_VERSION MATCHES "^OmpCluster")
        message(STATUS
            "Using OMPC Clang version ${CMAKE_C_COMPILER_VERSION}!")
    else()
        message(FATAL_ERROR
            "It seems you are using vanilla Clang instead of OMPC Clang. "
            "Please check your installation and try again!")
    endif()
else()
    message(FATAL_ERROR
        "Not using OMPC Clang as compiler, found: ${CMAKE_C_COMPILER_ID}.\n"
        ""
        "Make sure to export the following variables and call CMake from within"
        " the container.\n"
        ""
        "    export CC=clang\n"
        "    export CXX=clang++")
endif()

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# RelWithDebInfo use O3 by default to allow vectorization profiling
#option(RELWITHDEBINFO_FORCE_O3 "Force O3 with RelWithDebInfo" ON)
#if(RELWITHDEBINFO_FORCE_O3)
#  string(REPLACE "-O2" "-O3" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
#  string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
#endif()

# RelWithDebInfo use O3 by default to allow vectorization profiling
option(VERBOSE_VECTORIZE "Show vectorization report" OFF)
if(VERBOSE_VECTORIZE)
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")
  elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopt-info-vec-missed")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-missed")
  elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -opt-report5")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -opt-report5")
  else()
    message("Disable vectorization report: Unknown compiler ${CMAKE_C_COMPILER_ID}")
  endif()
endif()

option(USE_FASTMATH "Enable fast-math optimizations" ON)
if(NOT USE_FASTMATH)
  message("Manually disabled fast math compilation optimizations")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
else()
  message("Not building with fast-math optimizations: compiler not supported.")
endif()

option(USE_NATIVE "Build with -march=native" ON)
if(USE_NATIVE)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    message("Building with -march=native")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  else()
    message("Not building with -march=native: not supported by compiler")
  endif()
else()
  message("Not building with -march=native: not enabled.")
endif()

# Target: OMPCRuntime
# ------------------------------------------------------------------------------

add_library(OMPC::Runtime IMPORTED INTERFACE)

set_property(
    TARGET OMPC::Runtime
    PROPERTY
        INTERFACE_COMPILE_OPTIONS -fopenmp -fopenmp-targets=x86_64-pc-linux-gnu)

set_property(
    TARGET OMPC::Runtime
    PROPERTY
        INTERFACE_LINK_OPTIONS -fopenmp -fopenmp-targets=x86_64-pc-linux-gnu)

file(GLOB
  CHECK_SOURCE_FILES
  src/*.cpp
  src/headers/*.h
)

include_directories(src/headers)

add_executable(${PROJECT_NAME}
  src/main.cpp
  src/System.cpp
  src/Input.cpp
  src/RandomNumber.cpp
  src/Math.cpp
  src/Time.cpp
  src/MoleculeTarget.cpp
  src/GasBuffer.cpp
  src/Equipotential.cpp
  src/LinkedCell.cpp
  src/Force.cpp
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE OMPC::Runtime)
