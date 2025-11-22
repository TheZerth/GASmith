#include <benchmark/benchmark.h>

#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/ops/geometric.h"
#include "ga/versor.h"
#include "ga/rotor.h"

using namespace ga;
using namespace ga::ops;

static Multivector basisVec(const Algebra& alg, int axis) {
    Multivector v(alg);
    v.setComponent(Blade::getBasis(axis), 1.0f);
    return v;
}

// -----------------------------------------------------------------------------
// Versor benchmarks
// -----------------------------------------------------------------------------

static void BM_VersorApply_E3(benchmark::State& state) {
    Signature sig(3,0,0,true);
    Algebra alg(sig);

    Multivector e1 = basisVec(alg, 0);
    Multivector e2 = basisVec(alg, 1);
    Multivector v  = basisVec(alg, 2);

    Versor V(alg, geometricProduct(e2, e1));

    for (auto _ : state) {
        benchmark::DoNotOptimize(V.apply(v));
    }
}
BENCHMARK(BM_VersorApply_E3);

static void BM_VersorInverse_E3(benchmark::State& state) {
    Signature sig(3,0,0,true);
    Algebra alg(sig);

    Multivector e1 = basisVec(alg, 0);
    Multivector e2 = basisVec(alg, 1);

    Versor V(alg, geometricProduct(e2, e1));

    for (auto _ : state) {
        benchmark::DoNotOptimize(V.inverse());
    }
}
BENCHMARK(BM_VersorInverse_E3);

// -----------------------------------------------------------------------------
// Rotor benchmarks
// -----------------------------------------------------------------------------

static void BM_RotorApply_E3(benchmark::State& state) {
    Signature sig(3,0,0,true);
    Algebra alg(sig);

    Multivector e1 = basisVec(alg, 0);
    Multivector e2 = basisVec(alg, 1);

    Rotor R = Rotor::fromPlaneAngle(e1, e2, M_PI / 3.0f);

    Multivector v = basisVec(alg, 2);

    for (auto _ : state) {
        benchmark::DoNotOptimize(R.apply(v));
    }
}
BENCHMARK(BM_RotorApply_E3);

static void BM_RotorNormalize_E3(benchmark::State& state) {
    Signature sig(3,0,0,true);
    Algebra alg(sig);

    Multivector e1 = basisVec(alg, 0);
    Multivector e2 = basisVec(alg, 1);

    for (auto _ : state) {
        Rotor R = Rotor::fromPlaneAngle(e1, e2, M_PI / 4.0f);
        benchmark::DoNotOptimize(R);
    }
}
BENCHMARK(BM_RotorNormalize_E3);

// -----------------------------------------------------------------------------
// Optional: STA performance comparison
// -----------------------------------------------------------------------------

static void BM_RotorApply_STA(benchmark::State& state) {
    // STA signature (1,3)
    Signature sig(1,3,0,true);
    Algebra alg(sig);

    Multivector e0 = basisVec(alg, 0);
    Multivector e1 = basisVec(alg, 1);

    Rotor R = Rotor::fromPlaneAngle(e0, e1, 0.25f);

    Multivector v = basisVec(alg, 2);

    for (auto _ : state) {
        benchmark::DoNotOptimize(R.apply(v));
    }
}
BENCHMARK(BM_RotorApply_STA);

// End benchmark file
