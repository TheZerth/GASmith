#include <gtest/gtest.h>
#include <cmath>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/storageDense.h"
#include "ga/ops/blade.h"
#include "ga/ops/geometric.h"
#include "ga/ops/wedge.h"

// Types from your library
using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

using ga::ops::geometricProduct;
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
    // (1,3,0) → diag(+1,-1,-1,-1)
    return Signature(/*p=*/1, /*q=*/3, /*r=*/0, /*isRightHanded=*/true);
}

static inline Signature make_pga3d_sig() {
    // (+,+,+,0)
    return Signature(/*p=*/3, /*q=*/0, /*r=*/1, /*isRightHanded=*/true);
}

// Construct a pure basis vector multivector: scale * e_i
static Multivector make_vector(const Algebra& alg, int axisIndex, double scale) {
    Multivector mv(alg);
    BladeMask mask = Blade::getBasis(axisIndex);
    mv.setComponent(mask, static_cast<float>(scale));
    return mv;
}

// Construct a pure scalar multivector: s * 1
static Multivector make_scalar(const Algebra& alg, double s) {
    Multivector mv(alg);
    mv.setComponent(static_cast<BladeMask>(0), static_cast<float>(s)); // scalar mask = 0
    return mv;
}

static bool multivectorAlmostEqual(const Multivector& A,
                                   const Multivector& B,
                                   double eps = 1e-6)
{
    if (A.alg->dimensions != B.alg->dimensions) return false;

    const int dims = A.alg->dimensions;
    const std::size_t n = (1u << dims);
    for (std::size_t i = 0; i < n; ++i) {
        double da = A.storage[i];
        double db = B.storage[i];
        if (std::fabs(da - db) > eps) return false;
    }
    return true;
}

static void expectMultivectorAlmostEqual(const Multivector& A,
                                         const Multivector& B,
                                         double eps = 1e-6)
{
    EXPECT_TRUE(multivectorAlmostEqual(A, B, eps));
}

static void expectMultivectorAlmostEqualScaled(const Multivector& A,
                                               const Multivector& B,
                                               double scale,
                                               double eps = 1e-6)
{
    ASSERT_EQ(A.alg->dimensions, B.alg->dimensions);
    const int dims = A.alg->dimensions;
    const std::size_t n = (1u << dims);
    for (std::size_t i = 0; i < n; ++i) {
        double da = A.storage[i];
        double db = B.storage[i];
        EXPECT_NEAR(da, scale * db, eps);
    }
}

// --------------------- Wedge tests (Euclidean3) -----------------------------

// e_i ^ e_i = 0 in any metric
TEST(WedgeProduct, VectorSelfWedgeZero_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    const int dims = alg.dimensions;
    for (int i = 0; i < dims; ++i) {
        Multivector ei = make_vector(alg, i, 1.0);

        Multivector w = wedge(ei, ei);

        const std::size_t n = (1u << dims);
        for (std::size_t k = 0; k < n; ++k) {
            EXPECT_NEAR(w.storage[k], 0.0, 1e-6);
        }
    }
}

// Antisymmetry on vectors: e_i ^ e_j = - e_j ^ e_i for i != j
TEST(WedgeProduct, Antisymmetry_OnVectors_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    const int dims = alg.dimensions;
    for (int i = 0; i < dims; ++i) {
        for (int j = 0; j < dims; ++j) {
            Multivector ei = make_vector(alg, i, 1.0);
            Multivector ej = make_vector(alg, j, 1.0);

            Multivector wij = wedge(ei, ej);
            Multivector wji = wedge(ej, ei);

            if (i == j) {
                // Already tested above, but let's assert zero again here
                const std::size_t n = (1u << dims);
                for (std::size_t k = 0; k < n; ++k) {
                    EXPECT_NEAR(wij.storage[k], 0.0, 1e-6);
                    EXPECT_NEAR(wji.storage[k], 0.0, 1e-6);
                }
            } else {
                // e_i ^ e_j = - e_j ^ e_i
                expectMultivectorAlmostEqualScaled(wij, wji, -1.0);
            }
        }
    }
}

// Concrete check: e1 ^ e2 = e12, and linearity on the right: e1 ^ (e2 + e3) = e12 + e13
TEST(WedgeProduct, ConcreteVectorBivectors_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector e1 = make_vector(alg, 0, 1.0);
    Multivector e2 = make_vector(alg, 1, 1.0);
    Multivector e3 = make_vector(alg, 2, 1.0);

    // e1 ^ e2
    Multivector e1we2 = wedge(e1, e2);
    BladeMask e12Mask = Blade::getBasis(0) | Blade::getBasis(1);

    // Only e12 component should be 1, everything else 0
    const int dims = alg.dimensions;
    const std::size_t n = (1u << dims);
    for (std::size_t k = 0; k < n; ++k) {
        double coeff = e1we2.storage[k];
        if (k == e12Mask) {
            EXPECT_NEAR(coeff, 1.0, 1e-6);
        } else {
            EXPECT_NEAR(coeff, 0.0, 1e-6);
        }
    }

    // e1 ^ (e2 + e3) = e12 + e13
    Multivector e2_plus_e3(alg);
    e2_plus_e3.setComponent(Blade::getBasis(1), 1.0f); // e2
    e2_plus_e3.setComponent(Blade::getBasis(2), 1.0f); // e3

    Multivector w = wedge(e1, e2_plus_e3);

    BladeMask e13Mask = Blade::getBasis(0) | Blade::getBasis(2);
    for (std::size_t k = 0; k < n; ++k) {
        double coeff = w.storage[k];
        if (k == e12Mask || k == e13Mask) {
            EXPECT_NEAR(coeff, 1.0, 1e-6);
        } else {
            EXPECT_NEAR(coeff, 0.0, 1e-6);
        }
    }
}

// Associativity on three independent vectors:
// (e1 ^ e2) ^ e3 = e1 ^ (e2 ^ e3) = e123
TEST(WedgeProduct, Associativity_OnVectors_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector e1 = make_vector(alg, 0, 1.0);
    Multivector e2 = make_vector(alg, 1, 1.0);
    Multivector e3 = make_vector(alg, 2, 1.0);

    Multivector e1we2 = wedge(e1, e2);
    Multivector e2we3 = wedge(e2, e3);

    Multivector left  = wedge(e1we2, e3);
    Multivector right = wedge(e1, e2we3);

    expectMultivectorAlmostEqual(left, right);

    // And check the result is exactly +e123
    BladeMask e123Mask = Blade::getBasis(0) | Blade::getBasis(1) | Blade::getBasis(2);
    const int dims = alg.dimensions;
    const std::size_t n = (1u << dims);
    for (std::size_t k = 0; k < n; ++k) {
        double coeff = left.storage[k];
        if (k == e123Mask) {
            EXPECT_NEAR(coeff, 1.0, 1e-6);
        } else {
            EXPECT_NEAR(coeff, 0.0, 1e-6);
        }
    }
}

// Scalars should behave like scalar multiplication: (s) ^ A = s A and A ^ (s) = s A
TEST(WedgeProduct, ScalarLinearity_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    double s = 2.5;
    Multivector scalar = make_scalar(alg, s);
    Multivector e1     = make_vector(alg, 0, 1.0);

    Multivector left  = wedge(scalar, e1);
    Multivector right = wedge(e1, scalar);

    // Expected: s * e1
    Multivector expected(alg);
    expected.setComponent(Blade::getBasis(0), static_cast<float>(s));

    expectMultivectorAlmostEqual(left, expected);
    expectMultivectorAlmostEqual(right, expected);
}

// --------------------- Metric-independence for orthogonal vectors -----------------------------
//
// For orthogonal basis vectors (no contraction), the wedge product should be
// metric-independent: e_i ^ e_j has the same blade mask in any diagonal signature.
TEST(WedgeProduct, OrthogonalVectors_MetricIndependent) {
    // STA
    {
        Signature sig = make_sta_sig();
        Algebra   alg(sig);

        Multivector e0 = make_vector(alg, 0, 1.0);
        Multivector e1 = make_vector(alg, 1, 1.0);

        Multivector w = wedge(e0, e1);
        BladeMask expectedMask = Blade::getBasis(0) | Blade::getBasis(1);

        const int dims = alg.dimensions;
        const std::size_t n = (1u << dims);
        for (std::size_t k = 0; k < n; ++k) {
            double coeff = w.storage[k];
            if (k == expectedMask) {
                EXPECT_NEAR(coeff, 1.0, 1e-6);
            } else {
                EXPECT_NEAR(coeff, 0.0, 1e-6);
            }
        }
    }

    // PGA3D (null axis involved)
    {
        Signature sig = make_pga3d_sig();
        Algebra   alg(sig);

        int nullAxisIndex = 3;
        Multivector e1     = make_vector(alg, 0, 1.0);
        Multivector eInf   = make_vector(alg, nullAxisIndex, 1.0);

        Multivector w = wedge(e1, eInf);
        BladeMask expectedMask = Blade::getBasis(0) | Blade::getBasis(nullAxisIndex);

        const int dims = alg.dimensions;
        const std::size_t n = (1u << dims);
        for (std::size_t k = 0; k < n; ++k) {
            double coeff = w.storage[k];
            if (k == expectedMask) {
                EXPECT_NEAR(coeff, 1.0, 1e-6);
            } else {
                EXPECT_NEAR(coeff, 0.0, 1e-6);
            }
        }
    }
}