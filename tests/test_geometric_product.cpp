// tests/test_geometric_product.cpp
#include <gtest/gtest.h>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/multivector.h"
#include "ga/ops/geometric.h"
#include "ga/ops/blade.h"
#include "ga/algebra.h"

using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Multivector;
using ga::Algebra;

// If Multivector / Algebra / DenseStorage live in a namespace, add using declarations here.
// For now I assume they are in the global namespace, as in your multivector.h snippet.

// Helper: make basis vector e_i
static Blade make_basis(int axisIndex) {
    return Blade{Blade::getBasis(axisIndex), +1};
}

// Helper: Euclidean3 signature (3,0,0)
static Signature make_euclidean3_sig() {
    return Signature(/*p=*/3, /*q=*/0, /*r=*/0, /*isRightHanded=*/true);
}

// ------------------------
// Basis sanity tests
// ------------------------

TEST(BasisTests, CombineBlade_WedgeOrientation) {
    Blade e1 = make_basis(0);
    Blade e2 = make_basis(1);

    Blade e12 = Blade::combineBlade(e1, e2);
    Blade e21 = Blade::combineBlade(e2, e1);

    // e1 ^ e2 = + e12
    EXPECT_EQ(e12.sign, +1);
    EXPECT_TRUE(Blade::hasAxis(e12.mask, 0));
    EXPECT_TRUE(Blade::hasAxis(e12.mask, 1));

    // e2 ^ e1 = - e12
    EXPECT_EQ(e21.sign, -1);
    EXPECT_EQ(e21.mask, e12.mask);
}

// ------------------------
// Blade-level geometric product
// ------------------------

TEST(GeometricProductBlade, Euclidean3_Simple) {
    Signature sig = make_euclidean3_sig();

    Blade e1 = make_basis(0);
    Blade e2 = make_basis(1);

    // e1*e1 = +1
    {
        Blade r = ga::ops::geometricProductBlade(e1, e1, sig);
        EXPECT_TRUE(Blade::isScalarBasis(r));
        EXPECT_EQ(r.sign, +1);
    }

    // e1*e2 = + e12
    {
        Blade r = ga::ops::geometricProductBlade(e1, e2, sig);
        EXPECT_FALSE(Blade::isScalarBasis(r));
        EXPECT_TRUE(Blade::hasAxis(r.mask, 0));
        EXPECT_TRUE(Blade::hasAxis(r.mask, 1));
        EXPECT_EQ(r.sign, +1);
    }

    // e2*e1 = - e12
    {
        Blade r = ga::ops::geometricProductBlade(e2, e1, sig);
        EXPECT_FALSE(Blade::isScalarBasis(r));
        EXPECT_TRUE(Blade::hasAxis(r.mask, 0));
        EXPECT_TRUE(Blade::hasAxis(r.mask, 1));
        EXPECT_EQ(r.sign, -1);
    }
}

// ------------------------
// Multivector geometric product
// ------------------------


// Helper: create a simple multivector in 3D: mv = 1 + e1 + 2 e2
static Multivector make_simple_mv(const Algebra& alg) {
    Multivector mv(alg);
    mv.setComponent(static_cast<BladeMask>(0), 1.0);                    // scalar
    mv.setComponent(Blade::getBasis(0), 1.0);                            // e1
    mv.setComponent(Blade::getBasis(1), 2.0);                            // e2
    return mv;
}

TEST(GeometricProduct, Multivector_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra alg(sig);

    Multivector A = make_simple_mv(alg);
    Multivector B = make_simple_mv(alg);

    Multivector C = ga::ops::geometricProduct(A, B);

    // A = 1 + e1 + 2 e2
    // B = 1 + e1 + 2 e2
    //
    // Compute expected scalar part:
    // scalar = 1*1 + (e1*e1) + (2e2*2e2)
    //        = 1 + 1 + 4 = 6 (since in Euclidean3, e1^2 = e2^2 = +1)
    EXPECT_NEAR(C.storage[0], 6.0, 1e-12);

    // Expected vector part:
    // e1: from 1*e1 + e1*1 + (2e2)*(e1) + (e1)*(2e2) components
    // The cross terms e2*e1 and e1*e2 become bivectors, not vectors,
    // so vector e1 coefficient should be 1 + 1 = 2.
    BladeMask e1Mask = Blade::getBasis(0);
    EXPECT_NEAR(C.storage[e1Mask], 2.0, 1e-12);

    // Same reasoning for e2:
    BladeMask e2Mask = Blade::getBasis(1);
    // e2 coefficient: 2 (from A) + 2 (from B) = 4
    EXPECT_NEAR(C.storage[e2Mask], 4.0, 1e-12);
}
