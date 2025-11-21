#pragma once
#include "blade.h"
#include "../multivector.h"
// Implements a full clifford product of two n-dimensional generalized multivectors.
// Multivecors must share a metric

using ga::BladeMask;

namespace ga::ops {

    // A runtime callback type: decides whether a given (gradeA, gradeB, gradeR) term from the geometric product should be kept.
    using GradeFilterFn = bool (*)(int gradeA, int gradeB, int gradeR);

    // Do full geometric product, keep only terms where keep(gradeA, gradeB, gradeR) == true.
    inline Multivector geometricProductFiltered(const Multivector& A,
                                           const Multivector& B,
                                           GradeFilterFn keep)
    {
        const Algebra* algA = A.alg;
        const Algebra* algB = B.alg;

        if (!algA || !algB || algA != algB) {
            throw std::invalid_argument(
                "ga::ops::geometricProductCore: Multivectors must share the same Algebra");
        }

        const Algebra* alg = algA;
        const int dims = alg->dimensions;
        const int bladeCount = 1 << dims;

        Multivector result(*alg);

        for (int i = 0; i < bladeCount; ++i) {
            BladeMask maskA = static_cast<BladeMask>(i);
            double coeffA = A.component(maskA);
            if (coeffA == 0.0)
                continue;

            int gradeA = 0;
            if (keep) {
                gradeA = ga::Blade::getGrade(maskA);
            }

            for (int j = 0; j < bladeCount; ++j) {
                BladeMask maskB = static_cast<BladeMask>(j);
                double coeffB = B.component(maskB);
                if (coeffB == 0.0)
                    continue;

                int gradeB = 0;
                if (keep) {
                    gradeB = ga::Blade::getGrade(maskB);
                }

                Blade a{maskA, +1};
                Blade b{maskB, +1};

                Blade gp = geometricProductBlade(a, b, alg->signature);
                if (Blade::isZero(gp))
                    continue;

                if (keep) {
                    int gradeR = ga::Blade::getGrade(gp.mask);
                    if (!keep(gradeA, gradeB, gradeR))
                        continue;
                }

                double contrib = coeffA * coeffB * static_cast<double>(gp.sign);
                if (contrib == 0.0)
                    continue;

                double prev = result.component(gp.mask);
                result.setComponent(gp.mask, prev + contrib);
            }
        }

        return result;
    }

    inline Multivector geometricProduct(const Multivector& A, const Multivector& B) {
        return geometricProductFiltered(A, B, nullptr);
    }

}

