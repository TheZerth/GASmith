#include <gtest/gtest.h>

#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/ops/geometric.h"
#include "ga/ops/involutions.h"
#include "ga/ops/wedge.h"
#include "ga/versor.h"
#include "ga/rotor.h"

using namespace ga;
using namespace ga::ops;

// Helper: construct unit basis vector e_i in given algebra
static Multivector basisVec(const Algebra& alg, int axis) {
    Multivector v(alg);
    v.setComponent(Blade::getBasis(axis), 1.0f);
    return v;
}

// -----------------------------------------------------------------------------
// Versor Tests
// -----------------------------------------------------------------------------

TEST(Versor, InverseSandwichIdentity) {
    Signature sig = Signature(3, 0, 0, true);
    Algebra alg(sig);

    Multivector a = basisVec(alg, 0); // e1
    Multivector b = basisVec(alg, 1); // e2

    // Create any non-degenerate versor V = a b
    Multivector Vmv = geometricProduct(a, b);
    Versor V(alg, Vmv);

    Multivector invV = V.inverse();

    // Check V * V^{-1} ≈ 1
    Multivector id = geometricProduct(Vmv, invV);
    EXPECT_NEAR(id.component(0), 1.0f, 1e-6);

    // All non-scalar parts should vanish
    const int dims = alg.dimensions;
    const std::size_t N = (1u << dims);
    for (std::size_t i = 1; i < N; ++i) {
        EXPECT_NEAR(id.component(static_cast<BladeMask>(i)), 0.0f, 1e-6);
    }
}

TEST(Versor, ApplyEqualsSandwich) {
    Signature sig = Signature(3, 0, 0, true);
    Algebra alg(sig);

    Multivector e1 = basisVec(alg, 0);
    Multivector e2 = basisVec(alg, 1);

    // V = e2 e1
    Multivector Vmv = geometricProduct(e2, e1);
    Versor V(alg, Vmv);

    // Apply to e1
    Multivector applied = V.apply(e1);

    // Compute manually: V e1 V^{-1}
    Multivector invV = V.inverse();
    Multivector manual = geometricProduct(geometricProduct(Vmv, e1), invV);

    const int dims = alg.dimensions;
    const std::size_t N = (1u << dims);
    for (std::size_t i = 0; i < N; ++i) {
        EXPECT_NEAR(applied.component(i), manual.component(i), 1e-6);
    }
}

// -----------------------------------------------------------------------------
// Rotor Tests
// -----------------------------------------------------------------------------

TEST(Rotor, Normalization) {
    Signature sig = Signature(3, 0, 0, true);
    Algebra alg(sig);

    Multivector e1 = basisVec(alg, 0);
    Multivector e2 = basisVec(alg, 1);

    Rotor R = Rotor::fromPlaneAngle(e1, e2, M_PI / 3.0f); // 60 degrees

    // Check R ~R ≈ 1
    Multivector rrev = reverse(R.value());
    Multivector n2 = geometricProduct(R.value(), rrev);

    EXPECT_NEAR(n2.component(0), 1.0f, 1e-6);

    const int dims = alg.dimensions;
    const std::size_t N = (1u << dims);
    for (std::size_t i = 1; i < N; ++i) {
        EXPECT_NEAR(n2.component(static_cast<BladeMask>(i)), 0.0f, 1e-6);
    }
}

TEST(Rotor, RotatesE1ToE2_90Deg) {
    Signature sig = Signature(3, 0, 0, true);
    Algebra alg(sig);

    Multivector e1 = basisVec(alg, 0);
    Multivector e2 = basisVec(alg, 1);

    Rotor R = Rotor::fromPlaneAngle(e1, e2, M_PI / 2.0f); // 90°

    Multivector rotated = R.apply(e1);

    // Expect rotated ≈ e2
    EXPECT_NEAR(rotated.component(Blade::getBasis(1)), 1.0f, 1e-6);

    // All other components ≈ 0
    EXPECT_NEAR(rotated.component(Blade::getBasis(0)), 0.0f, 1e-6);
    EXPECT_NEAR(rotated.component(Blade::getBasis(2)), 0.0f, 1e-6);
}

// End test file
