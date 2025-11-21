#include <gtest/gtest.h>
#include <cmath>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/storageDense.h"
#include "ga/ops/blade.h"
#include "ga/ops/dual.h"

using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

using ga::ops::dual;

// --------------------- Helpers -----------------------------

static inline Signature make_euclidean3_sig() {
    // (+,+,+)
    return Signature(/*p=*/3, /*q=*/0, /*r=*/0, /*isRightHanded=*/true);
}

static inline Signature make_sta_sig() {
    // (1,3,0) → (+,-,-,-)
    return Signature(/*p=*/1, /*q=*/3, /*r=*/0, /*isRightHanded=*/true);
}

// Simple builders
static Multivector make_scalar(const Algebra& alg, double s) {
    Multivector mv(alg);
    mv.setComponent(static_cast<BladeMask>(0), static_cast<float>(s)); // scalar mask = 0
    return mv;
}

static Multivector make_vector(const Algebra& alg, int axisIndex, double scale) {
    Multivector mv(alg);
    BladeMask mask = Blade::getBasis(axisIndex);
    mv.setComponent(mask, static_cast<float>(scale));
    return mv;
}

// A test multivector in 3D, using all grades 0..3
// mv = 1
//    + 2 e1 + 3 e2 + 4 e3
//    + 5 e12 + 6 e13 + 7 e23
//    + 8 e123
static Multivector make_test_mv_E3(const Algebra& alg) {
    Multivector mv(alg);

    // grade 0
    mv.setComponent(static_cast<BladeMask>(0), 1.0f);

    // grade 1
    BladeMask e1 = Blade::getBasis(0);
    BladeMask e2 = Blade::getBasis(1);
    BladeMask e3 = Blade::getBasis(2);
    mv.setComponent(e1, 2.0f);
    mv.setComponent(e2, 3.0f);
    mv.setComponent(e3, 4.0f);

    // grade 2
    BladeMask e12 = static_cast<BladeMask>(e1 | e2);
    BladeMask e13 = static_cast<BladeMask>(e1 | e3);
    BladeMask e23 = static_cast<BladeMask>(e2 | e3);
    mv.setComponent(e12, 5.0f);
    mv.setComponent(e13, 6.0f);
    mv.setComponent(e23, 7.0f);

    // grade 3
    BladeMask e123 = static_cast<BladeMask>(e1 | e2 | e3);
    mv.setComponent(e123, 8.0f);

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

static double coeff(const Multivector& mv, BladeMask m) {
    return mv.component(m);
}

// --------------------- Basis mapping tests (Euclidean3) -------------------
//
// Using standard orientation e1^e2^e3 and the implementation:
//
// I_mask  = 0b111
// dual(e_m) = sign * e_(I_mask ^ m), with sign from geometricProductBlade(e_m, e_(I^m))
//
// Expected mappings in 3D Euclidean:
//   dual(1)   =  e123
//   dual(e1)  =  e23
//   dual(e2)  = -e13
//   dual(e3)  =  e12
//   dual(e12) =  e3
//   dual(e13) = -e2
//   dual(e23) =  e1
//   dual(e123)=  1
TEST(Dual, BasisMapping_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    BladeMask scalar = static_cast<BladeMask>(0);
    BladeMask e1     = Blade::getBasis(0);
    BladeMask e2     = Blade::getBasis(1);
    BladeMask e3     = Blade::getBasis(2);
    BladeMask e12    = static_cast<BladeMask>(e1 | e2);
    BladeMask e13    = static_cast<BladeMask>(e1 | e3);
    BladeMask e23    = static_cast<BladeMask>(e2 | e3);
    BladeMask e123   = static_cast<BladeMask>(e1 | e2 | e3);

    // dual(1) = e123
    {
        Multivector one = make_scalar(alg, 1.0);
        Multivector d   = dual(one);

        EXPECT_NEAR(coeff(d, e123), 1.0, 1e-6);
        EXPECT_NEAR(coeff(d, scalar), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e1), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e2), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e3), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e12), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e13), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e23), 0.0, 1e-6);
    }

    // dual(e1) = e23
    {
        Multivector v = make_vector(alg, 0, 1.0);
        Multivector d = dual(v);

        EXPECT_NEAR(coeff(d, e23), 1.0, 1e-6);
        EXPECT_NEAR(coeff(d, e1), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e2), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e3), 0.0, 1e-6);
    }

    // dual(e2) = -e13
    {
        Multivector v = make_vector(alg, 1, 1.0);
        Multivector d = dual(v);

        EXPECT_NEAR(coeff(d, e13), -1.0, 1e-6);
        EXPECT_NEAR(coeff(d, e1), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e2), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e3), 0.0, 1e-6);
    }

    // dual(e3) = e12
    {
        Multivector v = make_vector(alg, 2, 1.0);
        Multivector d = dual(v);

        EXPECT_NEAR(coeff(d, e12), 1.0, 1e-6);
        EXPECT_NEAR(coeff(d, e1), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e2), 0.0, 1e-6);
        EXPECT_NEAR(coeff(d, e3), 0.0, 1e-6);
    }

    // dual(e12) = e3
    {
        Multivector mv(alg);
        mv.setComponent(e12, 1.0f);
        Multivector d = dual(mv);

        EXPECT_NEAR(coeff(d, e3), 1.0, 1e-6);
    }

    // dual(e13) = -e2
    {
        Multivector mv(alg);
        mv.setComponent(e13, 1.0f);
        Multivector d = dual(mv);

        EXPECT_NEAR(coeff(d, e2), -1.0, 1e-6);
    }

    // dual(e23) = e1
    {
        Multivector mv(alg);
        mv.setComponent(e23, 1.0f);
        Multivector d = dual(mv);

        EXPECT_NEAR(coeff(d, e1), 1.0, 1e-6);
    }

    // dual(e123) = 1
    {
        Multivector mv(alg);
        mv.setComponent(e123, 1.0f);
        Multivector d = dual(mv);

        EXPECT_NEAR(coeff(d, scalar), 1.0, 1e-6);
    }
}

// --------------------- Involution property in 3D Euclidean ---------------
//
// In 3D Euclidean, star^2 acts as identity on all grades:
//   star^2(A) = A
TEST(Dual, InvolutionProperty_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);

    Multivector d  = dual(mv);
    Multivector dd = dual(d);

    expectMultivectorAlmostEqual(mv, dd);
}

// --------------------- Linearity ------------------------------------------
//
// dual(A + B) = dual(A) + dual(B)
TEST(Dual, Linearity_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector A = make_test_mv_E3(alg);

    Multivector B(alg);
    // Choose a simple B: e1 + e12
    BladeMask e1  = Blade::getBasis(0);
    BladeMask e2  = Blade::getBasis(1);
    BladeMask e12 = static_cast<BladeMask>(e1 | e2);
    B.setComponent(e1, 1.0f);
    B.setComponent(e12, 1.0f);

    Multivector AplusB(alg);
    const int dims = alg.dimensions;
    const std::size_t n = (1u << dims);
    for (std::size_t i = 0; i < n; ++i) {
        BladeMask m = static_cast<BladeMask>(i);
        double c = A.component(m) + B.component(m);
        if (c != 0.0) {
            AplusB.setComponent(m, static_cast<float>(c));
        }
    }

    Multivector dA      = dual(A);
    Multivector dB      = dual(B);
    Multivector dAplusB = dual(AplusB);

    // Check dual(A + B) == dual(A) + dual(B)
    Multivector sum_dA_dB(alg);
    for (std::size_t i = 0; i < n; ++i) {
        BladeMask m = static_cast<BladeMask>(i);
        double c = dA.component(m) + dB.component(m);
        if (c != 0.0) {
            sum_dA_dB.setComponent(m, static_cast<float>(c));
        }
    }

    expectMultivectorAlmostEqual(dAplusB, sum_dA_dB);
}

// --------------------- Sanity check: no crash on 4D STA -------------------
//
// We won't assert a particular mapping here, just that dual() runs and
// preserves the total squared magnitude pattern (rough sanity).
TEST(Dual, STA_NoCrashBasicUse) {
    Signature sigSTA = make_sta_sig();
    Algebra   algSTA(sigSTA);

    // Build a test multivector using the first 3 spatial axes (0,1,2)
    Multivector mv = make_test_mv_E3(algSTA);

    Multivector d = dual(mv);
    Multivector dd = dual(d);

    // At least check dd is finite and not all zero if mv was non-zero
    const int dims = algSTA.dimensions;
    const std::size_t n = (1u << dims);

    double sumMv2 = 0.0;
    double sumDd2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double a = mv.storage[i];
        double b = dd.storage[i];
        sumMv2 += a*a;
        sumDd2 += b*b;
    }

    EXPECT_GT(sumMv2, 0.0);
    EXPECT_GT(sumDd2, 0.0);
}
