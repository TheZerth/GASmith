#include <benchmark/benchmark.h>
#include "ga/basis.h"

using ga::Blade;
using ga::BladeMask;

// --- Benchmark: getGrade ---
static void BM_getGrade(benchmark::State& state) {
    BladeMask masks[] = {
        0b00000000,
        0b00000001,
        0b00000101,
        0b11111111
    };
    for (auto _ : state) {
        for (auto m : masks) {
            benchmark::DoNotOptimize(Blade::getGrade(m));
        }
    }
}
BENCHMARK(BM_getGrade);

// --- Benchmark: hasAxis ---
static void BM_hasAxis(benchmark::State& state) {
    BladeMask m = 0b10101010;
    for (auto _ : state) {
        for (int i = 0; i < MAX_DIMENSIONS; ++i) {
            benchmark::DoNotOptimize(Blade::hasAxis(m, i));
        }
    }
}
BENCHMARK(BM_hasAxis);

// --- Benchmark: makeBlade ---
static void BM_makeBlade(benchmark::State& state) {
    int basis1[1] = {1};
    int basis2[2] = {1, 3};
    int basis3[3] = {3, 1, 2};

    for (auto _ : state) {
        benchmark::DoNotOptimize(Blade::makeBlade(basis1, 1));
        benchmark::DoNotOptimize(Blade::makeBlade(basis2, 2));
        benchmark::DoNotOptimize(Blade::makeBlade(basis3, 3));
    }
}
BENCHMARK(BM_makeBlade);

// --- Benchmark: combineBlade ---
static void BM_combineBlade(benchmark::State& state) {
    Blade e1 = Blade{Blade::getBasis(0), +1};
    Blade e2 = Blade{Blade::getBasis(1), +1};
    Blade e3 = Blade{Blade::getBasis(2), +1};

    for (auto _ : state) {
        benchmark::DoNotOptimize(Blade::combineBlade(e1, e2));
        benchmark::DoNotOptimize(Blade::combineBlade(e2, e3));
        benchmark::DoNotOptimize(Blade::combineBlade(e1, e3));
    }
}
BENCHMARK(BM_combineBlade);
