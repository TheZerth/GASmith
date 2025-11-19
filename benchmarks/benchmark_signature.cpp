#include <benchmark/benchmark.h>
#include "ga/signature.h"

using ga::Signature;
using ga::Metric;
using ga::Mask;

// Helpers to build some canonical test data -----------------------------

// Build a Metric for 3D Euclidean: (3,0,0), axisCount = 3
static Metric make_metric_euclidean3() {
    Metric m{};
    m[0] = 1;
    m[1] = 1;
    m[2] = 1;
    return m;
}

// Build a Metric for STA (1,3,0): time-like positive, 3 space-like negative, axisCount = 4
static Metric make_metric_sta_13() {
    Metric m{};
    m[0] =  1; // time
    m[1] = -1; // space x
    m[2] = -1; // space y
    m[3] = -1; // space z
    return m;
}

// Build masks for Euclidean3: pMask has e0,e1,e2; qMask and rMask empty
static void make_masks_euclidean3(Mask& pMask, Mask& qMask, Mask& rMask) {
    pMask.fill(false);
    qMask.fill(false);
    rMask.fill(false);

    pMask[0] = true;
    pMask[1] = true;
    pMask[2] = true;
}

// Build masks for STA (1,3,0): pMask has e0; qMask has e1,e2,e3; rMask empty
static void make_masks_sta_13(Mask& pMask, Mask& qMask, Mask& rMask) {
    pMask.fill(false);
    qMask.fill(false);
    rMask.fill(false);

    pMask[0] = true;      // time
    qMask[1] = true;      // space x
    qMask[2] = true;      // space y
    qMask[3] = true;      // space z
}

// ----------------------------------------------------------------------
// Benchmarks: construct Signature in different ways
// ----------------------------------------------------------------------

// 1) From (p,q,r,isRightHanded) for Euclidean3 (3,0,0)
static void BM_Signature_FromCounts_Euclidean3(benchmark::State& state) {
    for (auto _ : state) {
        Signature sig(3, 0, 0, true);
        benchmark::DoNotOptimize(sig);
    }
}

// 2) From Metric + axisCount for Euclidean3
static void BM_Signature_FromMetric_Euclidean3(benchmark::State& state) {
    const Metric m = make_metric_euclidean3();

    for (auto _ : state) {
        Signature sig(m, 3, true);
        benchmark::DoNotOptimize(sig);
    }
}

// 3) From Masks for Euclidean3
static void BM_Signature_FromMasks_Euclidean3(benchmark::State& state) {
    Mask pMask{}, qMask{}, rMask{};
    make_masks_euclidean3(pMask, qMask, rMask);

    for (auto _ : state) {
        Signature sig(pMask, qMask, rMask, true);
        benchmark::DoNotOptimize(sig);
    }
}

// 4) From (p,q,r,isRightHanded) for STA (1,3,0)
static void BM_Signature_FromCounts_STA(benchmark::State& state) {
    for (auto _ : state) {
        Signature sig(1, 3, 0, true);
        benchmark::DoNotOptimize(sig);
    }
}

// 5) From Metric + axisCount for STA (1,3,0)
static void BM_Signature_FromMetric_STA(benchmark::State& state) {
    const Metric m = make_metric_sta_13();

    for (auto _ : state) {
        Signature sig(m, 4, true);
        benchmark::DoNotOptimize(sig);
    }
}

// 6) From Masks for STA (1,3,0)
static void BM_Signature_FromMasks_STA(benchmark::State& state) {
    Mask pMask{}, qMask{}, rMask{};
    make_masks_sta_13(pMask, qMask, rMask);

    for (auto _ : state) {
        Signature sig(pMask, qMask, rMask, true);
        benchmark::DoNotOptimize(sig);
    }
}

// ----------------------------------------------------------------------
// Register benchmarks
// ----------------------------------------------------------------------

BENCHMARK(BM_Signature_FromCounts_Euclidean3);
BENCHMARK(BM_Signature_FromMetric_Euclidean3);
BENCHMARK(BM_Signature_FromMasks_Euclidean3);

BENCHMARK(BM_Signature_FromCounts_STA);
BENCHMARK(BM_Signature_FromMetric_STA);
BENCHMARK(BM_Signature_FromMasks_STA);

// main() is provided by benchmark::benchmark_main via CMake link flags
