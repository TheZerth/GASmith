#include <benchmark/benchmark.h>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/multivector.h"
#include "ga/algebra.h"
#include "ga/ops/geometric.h"
#include "ga/ops/blade.h"

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

    // Fill a small number of coefficients so we get interesting products
    // Example:
    // mv = 1 + e1 + 2 e23
    mv.setComponent(0b0000, 1.0);      // scalar
    mv.setComponent(Blade::getBasis(0), 1.0); // e1
    mv.setComponent(Blade::getBasis(1), 2.0); // e2
    mv.setComponent(Blade::getBasis(2), 3.0); // e3
    
    // Add a bivector to ensure we hit higher-grade products
    BladeMask e23 = Blade::getBasis(1) | Blade::getBasis(2);
    mv.setComponent(e23, 2.5);

    return mv;
}

// ---------------------------------------------------------
// Benchmark: Multivector Geometric Product (Euclidean 3D)
// Signature: (3,0,0)
// ---------------------------------------------------------
static void BM_MV_geometric_Euclidean3(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, true);
    Algebra alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize( geometricProduct(A, B) );
    }
}
BENCHMARK(BM_MV_geometric_Euclidean3);

// ---------------------------------------------------------
// Benchmark: Multivector Geometric Product (STA)
// Signature: (1,3,0)
// ---------------------------------------------------------
static void BM_MV_geometric_STA(benchmark::State& state) {
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, true);
    Algebra alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize( geometricProduct(A, B) );
    }
}
BENCHMARK(BM_MV_geometric_STA);

// ---------------------------------------------------------
// Benchmark: Multivector Geometric Product (PGA 3D)
// Signature: (3,0,1)
// ---------------------------------------------------------
static void BM_MV_geometric_PGA3D(benchmark::State& state) {
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, true);
    Algebra alg{sig};

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    for (auto _ : state) {
        benchmark::DoNotOptimize( geometricProduct(A, B) );
    }
}
BENCHMARK(BM_MV_geometric_PGA3D);
