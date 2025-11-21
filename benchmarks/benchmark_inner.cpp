#include <benchmark/benchmark.h>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/multivector.h"
#include "ga/algebra.h"
#include "ga/ops/geometric.h"
#include "ga/ops/blade.h"
#include "ga/ops/inner.h"

using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

// ---------------------------------------------------------
// Helper: simple multivector with a few nonzero components
// ---------------------------------------------------------

static Multivector make_simple_mv(const Algebra& alg) {
    Multivector mv(alg);

    // mv = 1 + e1 + 2 e2 + 3 e3 + 2.5 e23  (when dims >= 3)
    mv.setComponent(static_cast<BladeMask>(0), 1.0f); // scalar

    if (alg.dimensions >= 1) {
        mv.setComponent(Blade::getBasis(0), 1.0f);     // e1
    }
    if (alg.dimensions >= 2) {
        mv.setComponent(Blade::getBasis(1), 2.0f);     // e2
    }
    if (alg.dimensions >= 3) {
        BladeMask b2 = Blade::getBasis(1);
        BladeMask b3 = Blade::getBasis(2);
        BladeMask e23 = static_cast<BladeMask>(b2 | b3);
        mv.setComponent(e23, 2.5f);                    // e23
    }

    return mv;
}

// ---------------------------------------------------------
// Euclidean3: Signature (3,0,0)
// ---------------------------------------------------------

static void BM_MV_inner_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::inner(A, B));
    }
}
BENCHMARK(BM_MV_inner_Euclidean3);

static void BM_MV_leftContraction_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::leftContraction(A, B));
    }
}
BENCHMARK(BM_MV_leftContraction_Euclidean3);

static void BM_MV_rightContraction_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::rightContraction(A, B));
    }
}
BENCHMARK(BM_MV_rightContraction_Euclidean3);

// ---------------------------------------------------------
// STA: Signature (1,3,0)
// ---------------------------------------------------------

static void BM_MV_inner_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::inner(A, B));
    }
}
BENCHMARK(BM_MV_inner_STA);

static void BM_MV_leftContraction_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::leftContraction(A, B));
    }
}
BENCHMARK(BM_MV_leftContraction_STA);

static void BM_MV_rightContraction_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::rightContraction(A, B));
    }
}
BENCHMARK(BM_MV_rightContraction_STA);

// ---------------------------------------------------------
// PGA3D: Signature (3,0,1)
// ---------------------------------------------------------

static void BM_MV_inner_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::inner(A, B));
    }
}
BENCHMARK(BM_MV_inner_PGA3D);

static void BM_MV_leftContraction_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::leftContraction(A, B));
    }
}
BENCHMARK(BM_MV_leftContraction_PGA3D);

static void BM_MV_rightContraction_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::rightContraction(A, B));
    }
}
BENCHMARK(BM_MV_rightContraction_PGA3D);
