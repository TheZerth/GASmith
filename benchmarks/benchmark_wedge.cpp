#include <benchmark/benchmark.h>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/ops/geometric.h"
#include "ga/ops/blade.h"
#include "ga/ops/wedge.h"

using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

// ---------------------------------------------------------
// Helper functions
// ---------------------------------------------------------

static inline Blade make_basis_vector(int i) {
    return Blade{Blade::getBasis(i), +1};
}

static inline Multivector make_simple_mv(const Algebra& alg) {
    Multivector mv(alg);

    // Same pattern as your geometric benchmark:
    // mv = 1 + e1 + 2 e2 + 3 e3 + 2.5 e23
    mv.setComponent(static_cast<BladeMask>(0b0000), 1.0f); // scalar
    mv.setComponent(Blade::getBasis(0), 1.0f);             // e1
    mv.setComponent(Blade::getBasis(1), 2.0f);             // e2
    mv.setComponent(Blade::getBasis(2), 3.0f);             // e3

    BladeMask e23 = static_cast<BladeMask>(
        Blade::getBasis(1) | Blade::getBasis(2)
    );
    mv.setComponent(e23, 2.5f);

    return mv;
}

// ---------------------------------------------------------
// Benchmark: Multivector Wedge Product (Euclidean 3D)
// Signature: (3,0,0)
// ---------------------------------------------------------
static void BM_MV_wedge_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::wedge(A, B));
    }
}
BENCHMARK(BM_MV_wedge_Euclidean3);

// ---------------------------------------------------------
// Benchmark: Multivector Wedge Product (STA)
// Signature: (1,3,0)
// ---------------------------------------------------------
static void BM_MV_wedge_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::wedge(A, B));
    }
}
BENCHMARK(BM_MV_wedge_STA);

// ---------------------------------------------------------
// Benchmark: Multivector Wedge Product (PGA 3D)
// Signature: (3,0,1)
// ---------------------------------------------------------
static void BM_MV_wedge_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra   alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::wedge(A, B));
    }
}
BENCHMARK(BM_MV_wedge_PGA3D);

