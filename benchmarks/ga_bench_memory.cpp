#include "ga_bench_memory.h"

#if defined(_WIN32)
  #include <windows.h>
  #include <psapi.h>
#elif defined(__APPLE__)
  #include <mach/mach.h>
#elif defined(__linux__)
  #include <sys/resource.h>
#else
  #include <cstdlib>
#endif

namespace ga_bench {

std::uint64_t current_rss_bytes() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    if (GetProcessMemoryInfo(GetCurrentProcess(),
                             reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc),
                             sizeof(pmc))) {
        return static_cast<std::uint64_t>(pmc.WorkingSetSize);
    }
    return 0;
#elif defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(),
                  MACH_TASK_BASIC_INFO,
                  reinterpret_cast<task_info_t>(&info),
                  &count) == KERN_SUCCESS) {
        return static_cast<std::uint64_t>(info.resident_size);
    }
    return 0;
#elif defined(__linux__)
    struct rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // ru_maxrss is kilobytes on Linux
        return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ull;
    }
    return 0;
#else
    return 0;
#endif
}

void GAMemoryManager::Start() {
    // Called at the beginning of each benchmark run
    start_rss_ = current_rss_bytes();
}

void GAMemoryManager::Stop(Result& result) {
    // Called at the end of each benchmark run
    const auto end_rss = current_rss_bytes();
    const std::int64_t delta =
        (end_rss >= start_rss_) ? static_cast<std::int64_t>(end_rss - start_rss_)
                                : 0;

    // Interpretation (document this for yourself):
    //   max_bytes_used        ≈ RSS at end of run
    //   total_allocated_bytes ≈ RSS delta (end - start)
    //   net_heap_growth       ≈ RSS delta
    result.num_allocs            = 0;                       // not tracked
    result.max_bytes_used        = static_cast<std::int64_t>(end_rss);
    result.total_allocated_bytes = delta;
    result.net_heap_growth       = delta;
    // result.memory_iterations is filled by Google Benchmark.
}

} // namespace ga_bench
