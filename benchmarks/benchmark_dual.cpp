#include <benchmark/benchmark.h>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/multivector.h"
#include "ga/algebra.h"
#include "ga/ops/dual.h"

using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;
using ga::ops::dual;

// ---------------------------------------------------------
// Helper: multivector with several non-zero grades
// ---------------------------------------------------------

static Multivector make_simple_mv(const Algebra& alg) {
    Multivector mv(alg);

    // scalar
    mv.setComponent(static_cast<BladeMask>(0), 1.0f);

    if (alg.dimensions >= 1) {
        mv.setComponent(Blade::getBasis(0), 2.0f); // e1
    }
    if (alg.dimensions >= 2) {
        mv.setComponent(Blade::getBasis(1), 3.0f); // e2
    }
    if (alg.dimensions >= 3) {
        mv.setComponent(Blade::getBasis(2), 4.0f); // e3

        BladeMask e12 = static_cast<BladeMask>(
            Blade::getBasis(0) | Blade::getBasis(1));
        BladeMask e13 = static_cast<BladeMask>(
            Blade::getBasis(0) | Blade::getBasis(2));
        BladeMask e23 = static_cast<BladeMask>(
            Blade::getBasis(1) | Blade::getBasis(2));

        mv.setComponent(e12, 5.0f);
        mv.setComponent(e13, 6.0f);
        mv.setComponent(e23, 7.0f);

        BladeMask e123 = static_cast<BladeMask>(
            Blade::getBasis(0) | Blade::getBasis(1) | Blade::getBasis(2));
        mv.setComponent(e123, 8.0f);
    }

    return mv;
}

// ---------------------------------------------------------
// Euclidean3: (3,0,0)
// ---------------------------------------------------------

static void BM_Dual_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(dual(A));
    }
}
BENCHMARK(BM_Dual_Euclidean3);

// ---------------------------------------------------------
// STA: (1,3,0)
// ---------------------------------------------------------

static void BM_Dual_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(dual(A));
    }
}
BENCHMARK(BM_Dual_STA);

// ---------------------------------------------------------
// PGA3D: (3,0,1)
// ---------------------------------------------------------

static void BM_Dual_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(dual(A));
    }
}
BENCHMARK(BM_Dual_PGA3D);
