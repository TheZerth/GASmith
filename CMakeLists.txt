cmake_minimum_required(VERSION 4.0)
project(GASmith LANGUAGES CXX)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------------------------
# GASmith static library
# ------------------------------------------------------------------------------

add_library(GASmith STATIC
        library.cpp
        include/ga/signature.h
        include/ga/policies.h
        include/ga/algebra.h
        include/ga/basis.h
        include/ga/storageDense.h
        include/ga/multivector.h
        include/ga/linearMap.h
        include/ga/versor.h
        include/ga/rotor.h
        include/ga/ops/geometric.h
        include/ga/ops/wedge.h
        include/ga/ops/inner.h
        include/ga/ops/involutions.h
        include/ga/ops/dual.h
)

# Public headers live in include/
target_include_directories(GASmith
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ------------------------------------------------------------------------------
# Google Benchmark (fetched via FetchContent)
# ------------------------------------------------------------------------------

include(FetchContent)

FetchContent_Declare(
        benchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG        main
)

# Disable benchmark's own tests
set(BENCHMARK_ENABLE_GTEST   OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(benchmark)

# ------------------------------------------------------------------------------
# Benchmark executable
# ------------------------------------------------------------------------------

# You will create this file: benchmarks/benchmark_signature.cpp
add_executable(GASmith_bench
        benchmarks/benchmark_signature.cpp
        benchmarks/benchmark_signature.cpp
        benchmarks/benchmark_signature.h
)

target_link_libraries(GASmith_bench
        PRIVATE
        GASmith
        benchmark::benchmark
        benchmark::benchmark_main
)

# ------------------------------------------------------------------------------
# Custom target to run benchmarks and store results
# ------------------------------------------------------------------------------

# Benchmarks will be written to: <build-dir>/benchmarks/latest.json
set(GASMITH_BENCH_OUTPUT_DIR ${CMAKE_BINARY_DIR}/benchmarks)
set(GASMITH_BENCH_OUTPUT_FILE ${GASMITH_BENCH_OUTPUT_DIR}/latest.json)

add_custom_target(run_benchmarks
        COMMAND ${CMAKE_COMMAND} -E make_directory ${GASMITH_BENCH_OUTPUT_DIR}
        COMMAND GASmith_bench --benchmark_format=json --benchmark_out=${GASMITH_BENCH_OUTPUT_FILE}
        DEPENDS GASmith_bench
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running GASmith benchmarks and writing results to ${GASMITH_BENCH_OUTPUT_FILE}"
)
