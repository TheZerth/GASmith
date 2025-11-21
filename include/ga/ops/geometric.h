#pragma once
#include "blade.h"
#include "../multivector.h"
// Implements a full clifford product of two n-dimensional generalized multivectors.
// Multivecors must share a metric

using ga::BladeMask;

namespace ga::ops {

    inline Multivector geometricProduct(const Multivector& A, const Multivector& B) {
    if (A.alg != B.alg) {
        throw std::invalid_argument("ga::ops::geometricProduct: Multivectors must share the same Algebra");
    }
    if (!A.alg || !B.alg) {
        throw std::invalid_argument("ga::ops::geometricProduct: Multivectors must have an Algebra");
    }
    Multivector out(*A.alg);
        for (BladeMask i = 0; i < (1u << A.alg->dimensions); ++i) {
            double a = A.storage[i];
            if (a == 0) continue;

            Blade ba{i, +1};

            for (BladeMask j = 0; j < (1u << A.alg->dimensions); ++j) {
                double b = B.storage[j];
                if (b == 0) continue;

                Blade bb{j, +1};
                Blade bc = ga::ops::geometricProductBlade(ba, bb, A.alg->signature);

                if (!Blade::isZero(bc)) {
                    out.storage[bc.mask] += a * b * bc.sign;
                }
            }
        }
        return out;
    }

    // A runtime callback type: decides whether a given (gradeA, gradeB, gradeR)
    // term from the geometric product should be kept.
    using GradeFilterFn = bool (*)(int gradeA, int gradeB, int gradeR);

    // Non-templated version: do full geometric product, keep only terms where
    // keep(gradeA, gradeB, gradeR) == true.
    inline Multivector geometricProductFiltered(const Multivector& A,
                                           const Multivector& B,
                                           GradeFilterFn keep)
    {
        const Algebra* algebraA = A.alg;
        const Algebra* algebraB = B.alg;

        if (!algebraA || !algebraB || algebraA != algebraB) {
            throw std::invalid_argument(
                "ga::ops::binaryGradeFiltered: Multivectors must share the same Algebra");
        }

        // Both algebras the same, copy for reference.
        const Algebra* resultAlgebra = algebraA;
        const int usedDimensions = resultAlgebra->dimensions;
        const int bladeCount = 1 << usedDimensions;

        Multivector resultMV(*resultAlgebra);

        for (int i = 0; i < bladeCount; ++i) {
            BladeMask maskA = static_cast<BladeMask>(i);
            double coeffA = A.component(maskA);
            if (coeffA == 0.0)
                continue;

            int gradeA = ga::Blade::getGrade(maskA);

            for (int j = 0; j < bladeCount; ++j) {
                BladeMask maskB = static_cast<BladeMask>(j);
                double coeffB = B.component(maskB);
                if (coeffB == 0.0)
                    continue;

                int gradeB = ga::Blade::getGrade(maskB);

                Blade bladeA{maskA, +1};
                Blade bladeB{maskB, +1};

                Blade resultBlade = geometricProductBlade(bladeA, bladeB, resultAlgebra->signature);
                int gradeR = ga::Blade::getGrade(resultBlade.mask);

                if (!keep(gradeA, gradeB, gradeR))
                    continue;

                double contrib = coeffA * coeffB * static_cast<double>(resultBlade.sign);
                if (contrib == 0.0)
                    continue;

                double prev = resultMV.component(resultBlade.mask);
                resultMV.setComponent(resultBlade.mask, prev + contrib);
            }
        }

        return resultMV;
    }

}

