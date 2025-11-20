#include <benchmark/benchmark.h>
#include <cstdlib>
#include <string>

#include "ga_bench_memory.h"

namespace {

    std::string getenv_or(const char* key, const char* fallback) {
        if (const char* v = std::getenv(key)) return v;
        return fallback;
    }

    struct GASmithBenchmarkGlobalInit {
        GASmithBenchmarkGlobalInit() {
            using benchmark::AddCustomContext;
            using benchmark::RegisterMemoryManager;

            // 1) Register MemoryManager
            static ga_bench::GAMemoryManager memoryManager;
            RegisterMemoryManager(&memoryManager);

            // 2) Custom context (will appear in JSON under "context")
            AddCustomContext("build_type",   getenv_or("GA_BENCH_BUILD_TYPE", "unknown"));
            AddCustomContext("compiler",     getenv_or("GA_BENCH_COMPILER", "unknown"));
            AddCustomContext("ga_signature", getenv_or("GA_BENCH_SIGNATURE", "unknown"));
            AddCustomContext("git_sha",      getenv_or("GA_BENCH_GIT_SHA", "unknown"));
            AddCustomContext("git_branch",   getenv_or("GA_BENCH_GIT_BRANCH", "unknown"));
            AddCustomContext("run_id",       getenv_or("GA_BENCH_RUN_ID", "unknown"));
        }
    };

    // Static instance: its ctor runs before benchmark_main's main()
    GASmithBenchmarkGlobalInit g_benchInit;

} // namespace
