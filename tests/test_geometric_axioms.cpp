#include <gtest/gtest.h>
#include <cmath>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/storageDense.h"
#include "ga/ops/blade.h"
#include "ga/ops/geometric.h"

// Types from your library
using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

// geometricProduct is a free function (global namespace) from geometric.h
using ga::ops::geometricProduct;

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

// Make a simple 3D multivector: A = a0 + a1 e1 + a2 e2
static Multivector make_mv_3d(const Algebra& alg,
                              double a0, double a1, double a2)
{
    Multivector mv(alg);
    // DenseStorage is zeroed in ctor, so we just set non-zero components
    mv.setComponent(static_cast<BladeMask>(0), static_cast<float>(a0));          // scalar
    mv.setComponent(Blade::getBasis(0), static_cast<float>(a1));                 // e1
    mv.setComponent(Blade::getBasis(1), static_cast<float>(a2));                 // e2
    return mv;
}

// --------------------- Blade-level axioms -----------------------------

// Clifford relation for diagonal metric:
// e_i e_j + e_j e_i = 2 g_ij.  For i != j and orthogonal basis → g_ij = 0,
// so e_i e_j = - e_j e_i.
TEST(BladeGeometricAxioms, CliffordRelation_OffDiagonal_Euclidean3) {
    Signature sig = make_euclidean3_sig();

    for (int i = 0; i < sig.dimensionsUsed(); ++i) {
        for (int j = 0; j < sig.dimensionsUsed(); ++j) {
            Blade ei = make_basis(i);
            Blade ej = make_basis(j);

            Blade gij1 = ga::ops::geometricProductBlade(ei, ej, sig);
            Blade gij2 = ga::ops::geometricProductBlade(ej, ei, sig);

            if (i == j) {
                // e_i^2 = g_ii = +1 in Euclidean3
                EXPECT_TRUE(Blade::isScalarBasis(gij1));
                EXPECT_EQ(gij1.sign, +1);
            } else {
                // Off-diagonal: e_i e_j = - e_j e_i
                EXPECT_EQ(gij1.mask, gij2.mask);
                EXPECT_EQ(gij1.sign, -gij2.sign);
            }
        }
    }
}

// STA diagonal metric signs: e0^2 = +1, e1^2 = e2^2 = e3^2 = -1
TEST(BladeGeometricAxioms, CliffordDiagonal_STA) {
    Signature sig = make_sta_sig();
    ASSERT_EQ(sig.dimensionsUsed(), 4);

    // e0 time-like
    {
        Blade e0 = make_basis(0);
        Blade r = ga::ops::geometricProductBlade(e0, e0, sig);
        EXPECT_TRUE(Blade::isScalarBasis(r));
        EXPECT_EQ(r.sign, +1);
    }

    // e1,e2,e3 spatial
    for (int i = 1; i < 4; ++i) {
        Blade ei = make_basis(i);
        Blade r = ga::ops::geometricProductBlade(ei, ei, sig);
        EXPECT_TRUE(Blade::isScalarBasis(r));
        EXPECT_EQ(r.sign, -1);
    }
}

// PGA3D: null axis squared = 0
TEST(BladeGeometricAxioms, NullAxis_PGA3D) {
    Signature sig = make_pga3d_sig();
    ASSERT_EQ(sig.dimensionsUsed(), 4);

    int nullAxisIndex = 3; // 3rd index is the null axis (metric 0)
    Blade eInf = make_basis(nullAxisIndex);
    Blade r = ga::ops::geometricProductBlade(eInf, eInf, sig);

    EXPECT_TRUE(Blade::isZero(r));
}

// Associativity on basis vectors: (e_i e_j) e_k = e_i (e_j e_k)
TEST(BladeGeometricAxioms, Associativity_OnVectors_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    const int dims = sig.dimensionsUsed();

    for (int i = 0; i < dims; ++i) {
        for (int j = 0; j < dims; ++j) {
            for (int k = 0; k < dims; ++k) {
                Blade ei = make_basis(i);
                Blade ej = make_basis(j);
                Blade ek = make_basis(k);

                Blade left  = ga::ops::geometricProductBlade(
                    ga::ops::geometricProductBlade(ei, ej, sig),
                    ek,
                    sig
                );

                Blade right = ga::ops::geometricProductBlade(
                    ei,
                    ga::ops::geometricProductBlade(ej, ek, sig),
                    sig
                );

                EXPECT_EQ(left.mask, right.mask);
                EXPECT_EQ(left.sign, right.sign);
            }
        }
    }
}

// --------------------- Multivector axioms -----------------------------

// Scalar identity: 1 is left/right identity for geometric product
TEST(MultivectorGeometricAxioms, ScalarIdentity_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra alg(sig);

    Multivector A = make_mv_3d(alg, 1.0, 2.0, -3.0);

    Multivector one(alg);
    // ctor zeroes storage; set scalar to 1
    one.setComponent(static_cast<BladeMask>(0), 1.0f);

    Multivector left  = geometricProduct(one, A);
    Multivector right = geometricProduct(A, one);

    expectMultivectorAlmostEqual(left, A);
    expectMultivectorAlmostEqual(right, A);
}

// Bilinearity: (A + B) C = AC + BC
TEST(MultivectorGeometricAxioms, Bilinearity_Left_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra alg(sig);

    Multivector A = make_mv_3d(alg, 1.0, 2.0, -1.0);
    Multivector B = make_mv_3d(alg, -0.5, 0.0, 4.0);
    Multivector C = make_mv_3d(alg, 0.25, -1.0, 3.0);

    const int dims = alg.dimensions;
    const std::size_t n = (1u << dims);

    // A + B
    Multivector AplusB(alg);
    for (std::size_t i = 0; i < n; ++i) {
        AplusB.storage[i] = A.storage[i] + B.storage[i];
    }

    Multivector lhs = geometricProduct(AplusB, C);
    Multivector rhs = geometricProduct(A, C);
    Multivector tmp = geometricProduct(B, C);
    for (std::size_t i = 0; i < n; ++i) {
        rhs.storage[i] += tmp.storage[i];
    }

    expectMultivectorAlmostEqual(lhs, rhs);
}

// Bilinearity: A (B + C) = AB + AC
TEST(MultivectorGeometricAxioms, Bilinearity_Right_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra alg(sig);

    Multivector A = make_mv_3d(alg, 1.0, 2.0, -1.0);
    Multivector B = make_mv_3d(alg, -0.5, 0.0, 4.0);
    Multivector C = make_mv_3d(alg, 0.25, -1.0, 3.0);

    const int dims = alg.dimensions;
    const std::size_t n = (1u << dims);

    // B + C
    Multivector BplusC(alg);
    for (std::size_t i = 0; i < n; ++i) {
        BplusC.storage[i] = B.storage[i] + C.storage[i];
    }

    Multivector lhs = geometricProduct(A, BplusC);
    Multivector rhs = geometricProduct(A, B);
    Multivector tmp = geometricProduct(A, C);
    for (std::size_t i = 0; i < n; ++i) {
        rhs.storage[i] += tmp.storage[i];
    }

    expectMultivectorAlmostEqual(lhs, rhs);
}

// Associativity: (AB)C = A(BC) for small multivectors
TEST(MultivectorGeometricAxioms, Associativity_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra alg(sig);

    Multivector A = make_mv_3d(alg, 1.0, 2.0, -1.0);
    Multivector B = make_mv_3d(alg, -0.5, 0.0, 4.0);
    Multivector C = make_mv_3d(alg, 0.25, -1.0, 3.0);

    Multivector AB = geometricProduct(A, B);
    Multivector BC = geometricProduct(B, C);

    Multivector left  = geometricProduct(AB, C);
    Multivector right = geometricProduct(A, BC);

    expectMultivectorAlmostEqual(left, right);
}

// Consistency: if A = e_i, B = e_j, MV product should match blade-level GP
TEST(MultivectorGeometricAxioms, Consistency_WithBladeGP_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra alg(sig);

    const int dims = alg.dimensions;
    const std::size_t n = (1u << dims);

    for (int i = 0; i < dims; ++i) {
        for (int j = 0; j < dims; ++j) {
            Blade ei = make_basis(i);
            Blade ej = make_basis(j);

            // Multivectors A = e_i, B = e_j
            Multivector A(alg);
            Multivector B(alg);
            // ctor zeroes storage; just set these components
            A.setComponent(ei.mask, static_cast<float>(ei.sign));
            B.setComponent(ej.mask, static_cast<float>(ej.sign));

            Multivector C = geometricProduct(A, B);
            Blade bProd = ga::ops::geometricProductBlade(ei, ej, sig);

            for (std::size_t k = 0; k < n; ++k) {
                double coeff = C.storage[k];
                if (Blade::isZero(bProd)) {
                    EXPECT_NEAR(coeff, 0.0, 1e-6);
                } else {
                    if (k == bProd.mask) {
                        EXPECT_NEAR(coeff, static_cast<double>(bProd.sign), 1e-6);
                    } else {
                        EXPECT_NEAR(coeff, 0.0, 1e-6);
                    }
                }
            }
        }
    }
}
