#include <gtest/gtest.h>
#include <cmath>

#include "ga/basis.h"
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/storageDense.h"
#include "ga/ops/blade.h"
#include "ga/ops/geometric.h"
#include "ga/ops/involutions.h"

using ga::Blade;
using ga::BladeMask;
using ga::Signature;
using ga::Algebra;
using ga::Multivector;

using ga::ops::reverse;
using ga::ops::gradeInvolution;
using ga::ops::cliffordConjugate;

// --------------------- Helpers -----------------------------

static inline Signature make_euclidean3_sig() {
    // (+,+,+)
    return Signature(/*p=*/3, /*q=*/0, /*r=*/0, /*isRightHanded=*/true);
}

static inline Signature make_sta_sig() {
    // (1,3,0) → (+,-,-,-)
    return Signature(/*p=*/1, /*q=*/3, /*r=*/0, /*isRightHanded=*/true);
}

// Construct a scalar + vector + bivector + trivector multivector in 3D
// mv = a0
//    + a1 e1 + a2 e2 + a3 e3
//    + a12 e12 + a13 e13 + a23 e23
//    + a123 e123
static Multivector make_test_mv_E3(const Algebra& alg) {
    Multivector mv(alg);

    // scalar
    mv.setComponent(static_cast<BladeMask>(0), 1.0f);      // a0 = 1

    // grade 1
    mv.setComponent(Blade::getBasis(0), 2.0f);             // a1 = 2
    mv.setComponent(Blade::getBasis(1), 3.0f);             // a2 = 3
    mv.setComponent(Blade::getBasis(2), 4.0f);             // a3 = 4

    // grade 2
    BladeMask e12 = static_cast<BladeMask>(Blade::getBasis(0) | Blade::getBasis(1));
    BladeMask e13 = static_cast<BladeMask>(Blade::getBasis(0) | Blade::getBasis(2));
    BladeMask e23 = static_cast<BladeMask>(Blade::getBasis(1) | Blade::getBasis(2));

    mv.setComponent(e12, 5.0f);                            // a12 = 5
    mv.setComponent(e13, 6.0f);                            // a13 = 6
    mv.setComponent(e23, 7.0f);                            // a23 = 7

    // grade 3 (pseudoscalar e123)
    BladeMask e123 = static_cast<BladeMask>(
            Blade::getBasis(0)
          | Blade::getBasis(1)
          | Blade::getBasis(2));
    mv.setComponent(e123, 8.0f);                           // a123 = 8

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

// Convenience: read a component as double
static double coeff(const Multivector& mv, BladeMask m) {
    return mv.component(m);
}

// --------------------- Reverse tests -----------------------------

// Check sign pattern for reverse in 3D Euclidean:
// grade 0: +
// grade 1: +
// grade 2: -
// grade 3: -
TEST(Involutions, ReverseSignPattern_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);
    Multivector rv = reverse(mv);

    BladeMask scalar = static_cast<BladeMask>(0);
    BladeMask e1 = Blade::getBasis(0);
    BladeMask e2 = Blade::getBasis(1);
    BladeMask e3 = Blade::getBasis(2);
    BladeMask e12 = static_cast<BladeMask>(e1 | e2);
    BladeMask e13 = static_cast<BladeMask>(e1 | e3);
    BladeMask e23 = static_cast<BladeMask>(e2 | e3);
    BladeMask e123 = static_cast<BladeMask>(e1 | e2 | e3);

    // Original coefficients
    double a0   = coeff(mv, scalar);
    double a1   = coeff(mv, e1);
    double a2   = coeff(mv, e2);
    double a3   = coeff(mv, e3);
    double a12  = coeff(mv, e12);
    double a13  = coeff(mv, e13);
    double a23  = coeff(mv, e23);
    double a123 = coeff(mv, e123);

    // Check signs
    EXPECT_NEAR(coeff(rv, scalar), a0, 1e-6);   // grade 0: +
    EXPECT_NEAR(coeff(rv, e1),     a1, 1e-6);   // grade 1: +
    EXPECT_NEAR(coeff(rv, e2),     a2, 1e-6);
    EXPECT_NEAR(coeff(rv, e3),     a3, 1e-6);

    EXPECT_NEAR(coeff(rv, e12),  -a12, 1e-6);   // grade 2: -
    EXPECT_NEAR(coeff(rv, e13),  -a13, 1e-6);
    EXPECT_NEAR(coeff(rv, e23),  -a23, 1e-6);

    EXPECT_NEAR(coeff(rv, e123), -a123, 1e-6);  // grade 3: -
}

// Reverse applied twice should be identity
TEST(Involutions, ReverseInvolutionProperty) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);

    Multivector rr = reverse(reverse(mv));

    expectMultivectorAlmostEqual(mv, rr);
}

// --------------------- Grade involution tests -----------------------------

// Grade involution sign pattern in 3D:
// grade 0: +
// grade 1: -
// grade 2: +
// grade 3: -
TEST(Involutions, GradeInvolutionSignPattern_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);
    Multivector gv = gradeInvolution(mv);

    BladeMask scalar = static_cast<BladeMask>(0);
    BladeMask e1 = Blade::getBasis(0);
    BladeMask e2 = Blade::getBasis(1);
    BladeMask e3 = Blade::getBasis(2);
    BladeMask e12 = static_cast<BladeMask>(e1 | e2);
    BladeMask e13 = static_cast<BladeMask>(e1 | e3);
    BladeMask e23 = static_cast<BladeMask>(e2 | e3);
    BladeMask e123 = static_cast<BladeMask>(e1 | e2 | e3);

    double a0   = coeff(mv, scalar);
    double a1   = coeff(mv, e1);
    double a2   = coeff(mv, e2);
    double a3   = coeff(mv, e3);
    double a12  = coeff(mv, e12);
    double a13  = coeff(mv, e13);
    double a23  = coeff(mv, e23);
    double a123 = coeff(mv, e123);

    EXPECT_NEAR(coeff(gv, scalar), a0, 1e-6);    // grade 0: +
    EXPECT_NEAR(coeff(gv, e1),    -a1, 1e-6);    // grade 1: -
    EXPECT_NEAR(coeff(gv, e2),    -a2, 1e-6);
    EXPECT_NEAR(coeff(gv, e3),    -a3, 1e-6);

    EXPECT_NEAR(coeff(gv, e12),   a12, 1e-6);    // grade 2: +
    EXPECT_NEAR(coeff(gv, e13),   a13, 1e-6);
    EXPECT_NEAR(coeff(gv, e23),   a23, 1e-6);

    EXPECT_NEAR(coeff(gv, e123), -a123, 1e-6);   // grade 3: -
}

// Grade involution applied twice is identity
TEST(Involutions, GradeInvolutionInvolutionProperty) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);

    Multivector gg = gradeInvolution(gradeInvolution(mv));

    expectMultivectorAlmostEqual(mv, gg);
}

// --------------------- Clifford conjugation tests -------------------------

// Clifford conjugation sign pattern in 3D:
// sign = (-1)^{r(r+1)/2}
// r=0 → +
// r=1 → -
// r=2 → -
// r=3 → +
TEST(Involutions, CliffordConjugateSignPattern_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);
    Multivector cv = cliffordConjugate(mv);

    BladeMask scalar = static_cast<BladeMask>(0);
    BladeMask e1 = Blade::getBasis(0);
    BladeMask e2 = Blade::getBasis(1);
    BladeMask e3 = Blade::getBasis(2);
    BladeMask e12 = static_cast<BladeMask>(e1 | e2);
    BladeMask e13 = static_cast<BladeMask>(e1 | e3);
    BladeMask e23 = static_cast<BladeMask>(e2 | e3);
    BladeMask e123 = static_cast<BladeMask>(e1 | e2 | e3);

    double a0   = coeff(mv, scalar);
    double a1   = coeff(mv, e1);
    double a2   = coeff(mv, e2);
    double a3   = coeff(mv, e3);
    double a12  = coeff(mv, e12);
    double a13  = coeff(mv, e13);
    double a23  = coeff(mv, e23);
    double a123 = coeff(mv, e123);

    EXPECT_NEAR(coeff(cv, scalar), a0, 1e-6);    // grade 0: +
    EXPECT_NEAR(coeff(cv, e1),    -a1, 1e-6);    // grade 1: -
    EXPECT_NEAR(coeff(cv, e2),    -a2, 1e-6);
    EXPECT_NEAR(coeff(cv, e3),    -a3, 1e-6);

    EXPECT_NEAR(coeff(cv, e12),  -a12, 1e-6);    // grade 2: -
    EXPECT_NEAR(coeff(cv, e13),  -a13, 1e-6);
    EXPECT_NEAR(coeff(cv, e23),  -a23, 1e-6);

    EXPECT_NEAR(coeff(cv, e123),  a123, 1e-6);   // grade 3: +
}

// Clifford conjugation applied twice is identity
TEST(Involutions, CliffordConjugateInvolutionProperty) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);

    Multivector cc = cliffordConjugate(cliffordConjugate(mv));

    expectMultivectorAlmostEqual(mv, cc);
}

// --------------------- Relationship tests ---------------------------------
//
// Clifford conjugation = reverse ∘ gradeInvolution (and also gradeInvolution ∘ reverse)
TEST(Involutions, CompositionRelationships_Euclidean3) {
    Signature sig = make_euclidean3_sig();
    Algebra   alg(sig);

    Multivector mv = make_test_mv_E3(alg);

    Multivector rg = reverse(gradeInvolution(mv));
    Multivector gr = gradeInvolution(reverse(mv));
    Multivector c  = cliffordConjugate(mv);

    expectMultivectorAlmostEqual(rg, c);
    expectMultivectorAlmostEqual(gr, c);
}

// Metric independence: involutions do not depend on the metric signs
TEST(Involutions, MetricIndependence_STA) {
    // STA: (1,3,0), 4D space
    Signature sigSTA = make_sta_sig();
    Algebra   algSTA(sigSTA);

    // Build a test multivector using only the first 3 axes (0,1,2)
    Multivector mvSTA = make_test_mv_E3(algSTA);

    Multivector rvSTA = reverse(mvSTA);
    Multivector gvSTA = gradeInvolution(mvSTA);
    Multivector cvSTA = cliffordConjugate(mvSTA);

    // Euclidean3: (3,0,0), 3D space
    Signature sigE = make_euclidean3_sig();
    Algebra   algE(sigE);

    Multivector mvE  = make_test_mv_E3(algE);
    Multivector rvE  = reverse(mvE);
    Multivector gvE  = gradeInvolution(mvE);
    Multivector cvE  = cliffordConjugate(mvE);

    // Compare only the blades that use axes 0,1,2 (i.e. masks < 2^3).
    const int dimsE = algE.dimensions;             // 3
    const std::size_t n = (1u << dimsE);           // 8

    for (std::size_t i = 0; i < n; ++i) {
        BladeMask m = static_cast<BladeMask>(i);

        EXPECT_NEAR(rvSTA.component(m), rvE.component(m), 1e-6);
        EXPECT_NEAR(gvSTA.component(m), gvE.component(m), 1e-6);
        EXPECT_NEAR(cvSTA.component(m), cvE.component(m), 1e-6);
    }
}

