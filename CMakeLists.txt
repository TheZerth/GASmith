cmake_minimum_required(VERSION 3.28)
project(GASmith LANGUAGES CXX)

# If user provided Python_EXECUTABLE, map it to Python3_EXECUTABLE for FindPython3
if(DEFINED Python_EXECUTABLE AND NOT DEFINED Python3_EXECUTABLE)
    set(Python3_EXECUTABLE "${Python_EXECUTABLE}")
endif()

# Ensure Python support for automated benchmarking
find_package(Python3 COMPONENTS Interpreter REQUIRED)
# Setup venv
set(GA_VENV_DIR "${CMAKE_BINARY_DIR}/.venv")

if (WIN32)
    # Windows puts executables in Scripts/
    set(GA_VENV_PYTHON "${GA_VENV_DIR}/Scripts/python.exe")
    set(GA_VENV_PIP    "${GA_VENV_DIR}/Scripts/pip.exe")
else()
    # Linux/macOS puts executables in bin/
    set(GA_VENV_PYTHON "${GA_VENV_DIR}/bin/python")
    set(GA_VENV_PIP    "${GA_VENV_DIR}/bin/pip")
endif()
# add modules to venv
add_custom_target(setup_venv
        # Step A: Create venv if it doesn't exist
        COMMAND ${CMAKE_COMMAND} -E echo "Checking Python Virtual Environment..."
        COMMAND ${Python3_EXECUTABLE} -m venv ${GA_VENV_DIR}

        # Step B: Install requirements (Use 'pip install -r requirements.txt' for more packages)
        COMMAND ${CMAKE_COMMAND} -E echo "Installing dependencies into .venv..."
        COMMAND ${GA_VENV_PYTHON} -m pip install --upgrade pip
        COMMAND ${GA_VENV_PYTHON} -m pip install influxdb3-python
        COMMAND ${GA_VENV_PYTHON} -m pip install pandas
        COMMAND ${GA_VENV_PYTHON} -m pip install python-dotenv
        COMMAND ${GA_VENV_PYTHON} -m pip install certifi

        # Step C: Metadata
        COMMENT "Ensuring .venv exists and packages are installed"
        VERBATIM
)

# Ensure Git for benchmark info
find_package(Git QUIET)

if(Git_FOUND)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GA_GIT_SHA
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GA_GIT_BRANCH
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    set(GA_GIT_SHA "nogit")
    set(GA_GIT_BRANCH "nogit")
endif()

# Default GA signature label for benchmarks (you can override in cmake-gui / CLI)
set(GA_DEFAULT_SIGNATURE "unknown" CACHE STRING "Default GA signature name for benchmarks")

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
        include/ga/ops/blade.h
)

# Public headers live in include/
target_include_directories(GASmith
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Google Unit Tests
include(FetchContent)

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(GASmith_tests
        tests/test_geometric_product.cpp
        tests/test_geometric_product.cpp
        tests/test_geometric_axioms.cpp
)

target_link_libraries(GASmith_tests
        PRIVATE
        GASmith
        GTest::gtest_main
)

target_include_directories(GASmith_tests
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

include(GoogleTest)
gtest_discover_tests(GASmith_tests)

add_custom_target(run_tests
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS GASmith_tests
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running GASmith unit tests"
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
        benchmarks/ga_bench_memory.cpp
        benchmarks/ga_bench_init.cpp
        benchmarks/benchmark_basis.cpp
        benchmarks/benchmark_blade.cpp
        benchmarks/benchmark_multivector.cpp
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

# Get timestamp
string(TIMESTAMP TIME_TAG "%Y-%m-%d_%H_%M_%S")

set(GASMITH_BENCH_OUTPUT_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/results/)
set(GASMITH_BENCH_OUTPUT_FILE ${GASMITH_BENCH_OUTPUT_DIR}bench-${CMAKE_BUILD_TYPE}-${TIME_TAG}.json)

add_custom_target(run_benchmarks
        COMMAND ${CMAKE_COMMAND} -E make_directory ${GASMITH_BENCH_OUTPUT_DIR}

        # Run benchmarks with additional GA_BENCH_* context variables
        COMMAND ${CMAKE_COMMAND} -E env
        GA_BENCH_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        GA_BENCH_COMPILER=${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}
        GA_BENCH_GIT_SHA=${GA_GIT_SHA}
        GA_BENCH_GIT_BRANCH=${GA_GIT_BRANCH}
        GA_BENCH_SIGNATURE=${GA_DEFAULT_SIGNATURE}
        GA_BENCH_RUN_ID=bench-${CMAKE_BUILD_TYPE}-${TIME_TAG}
        $<TARGET_FILE:GASmith_bench>
        --benchmark_out=${GASMITH_BENCH_OUTPUT_FILE}
        --benchmark_out_format=json

        # Upload the JSON to InfluxDB
        # Ensure we use the Python3_EXECUTABLE found by find_package if Python_EXECUTABLE isn't set
        COMMAND ${GA_VENV_PYTHON}
        ${CMAKE_SOURCE_DIR}/tools/upload_benchmarks_to_influx.py
        --json ${GASMITH_BENCH_OUTPUT_FILE}

        DEPENDS GASmith_bench setup_venv
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running benchmarks with InfluxDB upload (using local .venv)"
)