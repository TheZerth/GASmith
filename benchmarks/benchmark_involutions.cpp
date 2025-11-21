#include <benchmark/benchmark.h>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/multivector.h"
#include "ga/algebra.h"
#include "ga/ops/involutions.h"

using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

using ga::ops::reverse;
using ga::ops::gradeInvolution;
using ga::ops::cliffordConjugate;

// ---------------------------------------------------------
// Helper: simple "dense" multivector with several grades
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

        BladeMask e12 = static_cast<BladeMask>(Blade::getBasis(0) | Blade::getBasis(1));
        BladeMask e13 = static_cast<BladeMask>(Blade::getBasis(0) | Blade::getBasis(2));
        BladeMask e23 = static_cast<BladeMask>(Blade::getBasis(1) | Blade::getBasis(2));
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
// Euclidean3: Signature (3,0,0)
// ---------------------------------------------------------

static void BM_Reverse_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(reverse(A));
    }
}
BENCHMARK(BM_Reverse_Euclidean3);

static void BM_GradeInvolution_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(gradeInvolution(A));
    }
}
BENCHMARK(BM_GradeInvolution_Euclidean3);

static void BM_CliffordConjugate_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(cliffordConjugate(A));
    }
}
BENCHMARK(BM_CliffordConjugate_Euclidean3);

// ---------------------------------------------------------
// STA: Signature (1,3,0)
// ---------------------------------------------------------

static void BM_Reverse_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(reverse(A));
    }
}
BENCHMARK(BM_Reverse_STA);

static void BM_GradeInvolution_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(gradeInvolution(A));
    }
}
BENCHMARK(BM_GradeInvolution_STA);

static void BM_CliffordConjugate_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(cliffordConjugate(A));
    }
}
BENCHMARK(BM_CliffordConjugate_STA);

// ---------------------------------------------------------
// PGA3D: Signature (3,0,1)
// ---------------------------------------------------------

static void BM_Reverse_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(reverse(A));
    }
}
BENCHMARK(BM_Reverse_PGA3D);

static void BM_GradeInvolution_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(gradeInvolution(A));
    }
}
BENCHMARK(BM_GradeInvolution_PGA3D);

static void BM_CliffordConjugate_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(cliffordConjugate(A));
    }
}
BENCHMARK(BM_CliffordConjugate_PGA3D);
