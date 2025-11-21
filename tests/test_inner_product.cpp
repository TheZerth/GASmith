#include <gtest/gtest.h>
#include <cmath>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/storageDense.h"
#include "ga/signature.h"
#include "ga/ops/blade.h"
#include "ga/ops/geometric.h"
#include "ga/ops/wedge.h"
#include "ga/ops/inner.h"

// Types from your library
using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

using ga::ops::inner;
using ga::ops::leftContraction;
using ga::ops::rightContraction;
using ga::ops::wedge;

// --------------------- Helpers -----------------------------

static inline Blade make_basis(int axisIndex) {
    return Blade{Blade::getBasis(axisIndex), +1};
}

static inline Signature make_euclidean3_sig() {
    // (+,+,+)
    return Signature(/*p=*/3, /*q=*/0, /*r=*/0, /*isRightHanded=*/true);
}

static inline Signature make_sta_sig() {
    // (1,3,0) → (+,-,-,-)
    return Signature(/*p=*/1, /*q=*/3, /*r=*/0, /*isRightHanded=*/true);
}

static inline Signature make_pga3d_sig() {
    // (+,+,+,0)
    return Signature(/*p=*/3, /*q=*/0, /*r=*/1, /*isRightHanded=*/true);
}

// Construct a pure basis-vector multivector: scale * e_i
static Multivector make_vector(const Algebra& alg, int axisIndex, double scale) {
    Multivector mv(alg);
    BladeMask m = Blade::getBasis(axisIndex);
    mv.setComponent(m, static_cast<float>(scale));
    return mv;
}

// Construct a scalar multivector: s * 1
static Multivector make_scalar(const Algebra& alg, double s) {
    Multivector mv(alg);
    mv.setComponent(static_cast<BladeMask>(0), static_cast<float>(s)); // scalar mask = 0
    return mv;
}

static bool multivectorAlmostEqual(const Multivector& A,
                                   const Multivector& B,
                                   double eps = 1e-6)
{
    if (!A.alg || !B.alg || A.alg != B.alg)
        return false;

    const int dims = A.alg->dimensions;
    const std::size_t n = (1u << dims);
    for (std::size_t i = 0; i < n; ++i) {
        double da = A.storage[i];
        double db = B.storage[i];
        if (std::fabs(da - db) > eps)
            return false;
    }
    return true;
}

static void expectMultivectorAlmostEqual(const Multivector& A,
                                         const Multivector& B,
                                         double eps = 1e-6)
{
    EXPECT_TRUE(multivectorAlmostEqual(A, B, eps));
}

static void expectScalarEqual(const Multivector& mv,
                              double expected,
                              double eps = 1e-6)
{
    ASSERT_TRUE(mv.alg != nullptr);
    const int dims = mv.alg->dimensions;
    const std::size_t n = (1u << dims);

    for (std::size_t k = 0; k < n; ++k) {
        double coeff = mv.storage[k];
        if (k == 0) { // scalar slot
            EXPECT_NEAR(coeff, expected, eps);
        } else {
            EXPECT_NEAR(coeff, 0.0, eps);
        }
    }
}

// --------------------- Inner product tests -----------------------------

// e_i ⋅ e_j = δ_ij in Euclidean3
TEST(InnerProduct, VectorDot_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    const int dims = alg.dimensions;
    for (int i = 0; i < dims; ++i) {
        for (int j = 0; j < dims; ++j) {
            Multivector ei = make_vector(alg, i, 1.0);
            Multivector ej = make_vector(alg, j, 1.0);

            Multivector dot = inner(ei, ej);

            double expected = (i == j) ? 1.0 : 0.0;
            expectScalarEqual(dot, expected);
        }
    }
}

// Inner product of vectors is symmetric in Euclidean signature
TEST(InnerProduct, Symmetry_Vectors_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    const int dims = alg.dimensions;
    for (int i = 0; i < dims; ++i) {
        for (int j = 0; j < dims; ++j) {
            Multivector ei = make_vector(alg, i, 1.0);
            Multivector ej = make_vector(alg, j, 1.0);

            Multivector dot_ij = inner(ei, ej);
            Multivector dot_ji = inner(ej, ei);

            expectMultivectorAlmostEqual(dot_ij, dot_ji);
        }
    }
}

// Check STA metric signs: e0⋅e0=+1, spatial axes negative
TEST(InnerProduct, VectorDot_STA) {
    Signature sig = make_sta_sig();
    Algebra   alg(sig);

    // e0 time-like, e1..e3 space-like under (1,3,0)
    for (int i = 0; i < 4; ++i) {
        Multivector ei = make_vector(alg, i, 1.0);
        Multivector dot = inner(ei, ei);

        double expected = (i == 0) ? +1.0 : -1.0;
        expectScalarEqual(dot, expected);
    }
}

// Inner product of a sum of vectors: (e1 + e2) ⋅ (e1 + e2) = 2
TEST(InnerProduct, SumOfVectors_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector e1 = make_vector(alg, 0, 1.0);
    Multivector e2 = make_vector(alg, 1, 1.0);

    Multivector a(alg);
    a.setComponent(Blade::getBasis(0), 1.0f); // e1
    a.setComponent(Blade::getBasis(1), 1.0f); // e2

    Multivector b = a;

    Multivector dot = inner(a, b);
    expectScalarEqual(dot, 2.0);
}

// --------------------- Contraction vs Inner on vectors --------------------
//
// For vectors, inner, leftContraction, and rightContraction agree:
//   e_i ⋅ e_j = e_i ⟍ e_j = e_i ⟎ e_j
TEST(Contraction, AgreesWithInner_OnVectors_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    const int dims = alg.dimensions;
    for (int i = 0; i < dims; ++i) {
        for (int j = 0; j < dims; ++j) {
            Multivector ei = make_vector(alg, i, 1.0);
            Multivector ej = make_vector(alg, j, 1.0);

            Multivector dot = inner(ei, ej);
            Multivector lc  = leftContraction(ei, ej);
            Multivector rc  = rightContraction(ei, ej);

            expectMultivectorAlmostEqual(dot, lc);
            expectMultivectorAlmostEqual(dot, rc);
        }
    }
}

// --------------------- Left contraction tests -----------------------------
//
// In Euclidean3:
//   e1 ⟍ (e1 ^ e2) = e2
//   e2 ⟍ (e1 ^ e2) = -e1
//   e3 ⟍ (e1 ^ e2) = 0
TEST(Contraction, LeftContraction_VectorBivector_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector e1 = make_vector(alg, 0, 1.0);
    Multivector e2 = make_vector(alg, 1, 1.0);
    Multivector e3 = make_vector(alg, 2, 1.0);

    Multivector biv = wedge(e1, e2); // e1 ^ e2

    Multivector c1 = leftContraction(e1, biv);
    Multivector c2 = leftContraction(e2, biv);
    Multivector c3 = leftContraction(e3, biv);

    // Expected: c1 = e2
    Multivector expected_c1(alg);
    expected_c1.setComponent(Blade::getBasis(1), 1.0f);

    // Expected: c2 = -e1
    Multivector expected_c2(alg);
    expected_c2.setComponent(Blade::getBasis(0), -1.0f);

    // Expected: c3 = 0
    Multivector expected_c3(alg); // all zero

    expectMultivectorAlmostEqual(c1, expected_c1);
    expectMultivectorAlmostEqual(c2, expected_c2);
    expectMultivectorAlmostEqual(c3, expected_c3);
}

// --------------------- Right contraction tests ----------------------------
//
// (e1 ^ e2) ⟎ e2 = e1
// (e1 ^ e2) ⟎ e1 = -e2
// (e1 ^ e2) ⟎ e3 = 0
TEST(Contraction, RightContraction_BivectorVector_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector e1 = make_vector(alg, 0, 1.0);
    Multivector e2 = make_vector(alg, 1, 1.0);
    Multivector e3 = make_vector(alg, 2, 1.0);

    Multivector biv = wedge(e1, e2); // e1 ^ e2

    Multivector c1 = rightContraction(biv, e2);
    Multivector c2 = rightContraction(biv, e1);
    Multivector c3 = rightContraction(biv, e3);

    // Expected: c1 = e1
    Multivector expected_c1(alg);
    expected_c1.setComponent(Blade::getBasis(0), 1.0f);

    // Expected: c2 = -e2
    Multivector expected_c2(alg);
    expected_c2.setComponent(Blade::getBasis(1), -1.0f);

    // Expected: c3 = 0
    Multivector expected_c3(alg); // all zero

    expectMultivectorAlmostEqual(c1, expected_c1);
    expectMultivectorAlmostEqual(c2, expected_c2);
    expectMultivectorAlmostEqual(c3, expected_c3);
}

// Grade-lowering behaviour: contracting a bivector with a scalar or
// scalar with a bivector should give 0 for left/right contraction.
TEST(Contraction, GradeLoweringBehaviour) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector scalar = make_scalar(alg, 3.14);
    Multivector e1 = make_vector(alg, 0, 1.0);
    Multivector e2 = make_vector(alg, 1, 1.0);

    Multivector biv = wedge(e1, e2);

    Multivector lc1 = leftContraction(biv, scalar);
    Multivector rc1 = rightContraction(scalar, biv);

    const int dims = alg.dimensions;
    const std::size_t n = (1u << dims);
    for (std::size_t k = 0; k < n; ++k) {
        EXPECT_NEAR(lc1.storage[k], 0.0, 1e-6);
        EXPECT_NEAR(rc1.storage[k], 0.0, 1e-6);
    }
}
