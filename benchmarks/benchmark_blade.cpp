#include <benchmark/benchmark.h>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/ops/blade.h"  // <-- declares ga::geometricProductBlade

using ga::Blade;
using ga::BladeMask;
using ga::Signature;

// Helper to make a basis vector e_k
static inline Blade make_basis_vector(int axisIndex) {
    return Blade{Blade::getBasis(axisIndex), +1};
}

// Helper to make a simple bivector e_i ^ e_j with i < j
static inline Blade make_bivector(int i, int j) {
    Blade ei = make_basis_vector(i);
    Blade ej = make_basis_vector(j);
    return Blade::combineBlade(ei, ej);
}

// Helper to make a simple trivector e_i ^ e_j ^ e_k with i < j < k
static inline Blade make_trivector(int i, int j, int k) {
    Blade ei = make_basis_vector(i);
    Blade ej = make_basis_vector(j);
    Blade ek = make_basis_vector(k);
    return Blade::combineBlade(Blade::combineBlade(ei, ej), ek);
}

// ---------------------------------------------------------
// Euclidean 3D: Signature (3, 0, 0)
// ---------------------------------------------------------
static void BM_geometricProductBlade_Euclidean3(benchmark::State& state) {
    // Right-handed 3D Euclidean: (+,+,+)
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/0, /*isRightHanded=*/true);

    Blade e1   = make_basis_vector(0);
    Blade e2   = make_basis_vector(1);
    Blade e3   = make_basis_vector(2);
    Blade e12  = make_bivector(0, 1);
    Blade e23  = make_bivector(1, 2);
    Blade e31  = make_bivector(2, 0);
    Blade e123 = make_trivector(0, 1, 2);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e1, e1, sig));   // e1*e1 = +1
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e1, e2, sig));   // e1*e2 = e12
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e2, e1, sig));   // e2*e1 = -e12
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e12, e3, sig));  // bivector * vector
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e3, e12, sig));  // vector * bivector
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e12, e23, sig)); // bivector * bivector
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e123, e1, sig)); // pseudoscalar * vector
    }
}
BENCHMARK(BM_geometricProductBlade_Euclidean3);

// ---------------------------------------------------------
// STA: Signature (1, 3, 0)  (1 time, 3 space)
// ---------------------------------------------------------
static void BM_geometricProductBlade_STA(benchmark::State& state) {
    // Common STA convention: (-,+,+,+) or (+,-,-,-).
    // Your Signature currently assigns +1 axes first, then -1 axes.
    // So (p=1, q=3, r=0) gives metric diag(+1, -1, -1, -1).
    Signature sig(/*p=*/1, /*q=*/3, /*r=*/0, /*isRightHanded=*/true);

    // Let e0 = time, e1,e2,e3 = space
    Blade e0   = make_basis_vector(0);
    Blade e1   = make_basis_vector(1);
    Blade e2   = make_basis_vector(2);
    Blade e3   = make_basis_vector(3);
    Blade e01  = make_bivector(0, 1);
    Blade e23  = make_bivector(2, 3);
    Blade e0123 = make_trivector(0, 1, 2); // not full 4D pseudoscalar, but still a good test

    for (auto _ : state) {
        // Time-like vs spatial-like products (mix of +1 and -1 metric entries)
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e0, e0, sig));   // e0*e0 = +1
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e1, e1, sig));   // e1*e1 = -1
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e0, e1, sig));   // e0*e1
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e1, e0, sig));   // e1*e0

        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e01, e23, sig)); // bivector * bivector
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e0123, e0, sig));// trivector * vector
    }
}
BENCHMARK(BM_geometricProductBlade_STA);

// ---------------------------------------------------------
// PGA 3D: Signature (3, 0, 1)  (3 Euclidean + 1 null)
// ---------------------------------------------------------
static void BM_geometricProductBlade_PGA3D(benchmark::State& state) {
    // PGA(3,0,1): (+,+,+,0)
    Signature sig(/*p=*/3, /*q=*/0, /*r=*/1, /*isRightHanded=*/true);

    // e0,e1,e2: Euclidean directions, e3: ideal / null direction
    Blade e0   = make_basis_vector(0);
    Blade e1   = make_basis_vector(1);
    Blade e2   = make_basis_vector(2);
    Blade eInf = make_basis_vector(3); // null basis

    Blade e01  = make_bivector(0, 1);
    Blade e2Inf = make_bivector(2, 3);

    for (auto _ : state) {
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e0, e0, sig));     // e0*e0 = +1
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(eInf, eInf, sig)); // eInf*eInf = 0
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e0, eInf, sig));   // mixed null / Euclidean
        benchmark::DoNotOptimize(ga::ops::geometricProductBlade(e01, e2Inf, sig)); // bivector * bivector
    }
}
BENCHMARK(BM_geometricProductBlade_PGA3D);

