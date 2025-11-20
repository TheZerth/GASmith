#pragma once

#include <cstdint>
#include <benchmark/benchmark.h>

namespace ga_bench {

    // Process-level RSS snapshot in bytes (approximate).
    std::uint64_t current_rss_bytes();

    // Memory manager that uses RSS snapshots per benchmark run.
    // It fills Google Benchmark's Result struct so JSON will contain:
    //   "max_bytes_used", "total_allocated_bytes", "net_heap_growth"
    class GAMemoryManager : public benchmark::MemoryManager {
    public:
        void Start() override;
        void Stop(Result& result) override;

    private:
        std::uint64_t start_rss_ = 0;
    };

} // namespace ga_bench
