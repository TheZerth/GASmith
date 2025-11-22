#pragma once
#include <stdexcept>

#include "blade.h"

// Implements a full clifford product of two n-dimensional generalized multivectors.
// Multivecors must share a metric

namespace ga::ops {

    using ga::BladeMask;
    using ga::Algebra;
    using ga::Multivector;
    using ga::Blade;

    // A runtime callback type: decides whether a given (gradeA, gradeB, gradeR) term from the geometric product should be kept.
    using GradeFilterFn = bool (*)(int gradeA, int gradeB, int gradeR);

    // Do full geometric product, keep only terms where keep(gradeA, gradeB, gradeR) == true.
    inline Multivector geometricProductFiltered(const Multivector& A,
                                           const Multivector& B,
                                           const GradeFilterFn keep)
    {
        const Algebra* algA = A.alg;

        if (const Algebra* algB = B.alg; !algA || !algB || algA != algB) {
            throw std::invalid_argument(
                "ga::ops::geometricProductFiltered: Multivectors must share the same Algebra");
        }

        const Algebra* alg = algA;
        const int dims = alg->dimensions;
        const int bladeCount = 1 << dims;

        Multivector result(*alg);

        for (int i = 0; i < bladeCount; ++i) {
            const auto maskA = static_cast<BladeMask>(i);
            const double coeffA = A.component(maskA);
            if (coeffA == 0.0)
                continue;

            int gradeA = 0;
            if (keep) {
                gradeA = ga::Blade::getGrade(maskA);
            }

            for (int j = 0; j < bladeCount; ++j) {
                const auto maskB = static_cast<BladeMask>(j);
                const double coeffB = B.component(maskB);
                if (coeffB == 0.0)
                    continue;

                int gradeB = 0;
                if (keep) {
                    gradeB = ga::Blade::getGrade(maskB);
                }

                Blade a{maskA, +1};
                Blade b{maskB, +1};

                const Blade gp = geometricProductBlade(a, b, alg->signature);
                if (Blade::isZero(gp))
                    continue;

                if (keep) {
                    if (const int gradeR = ga::Blade::getGrade(gp.mask); !keep(gradeA, gradeB, gradeR))
                        continue;
                }

                const double contrib = coeffA * coeffB * static_cast<double>(gp.sign);
                if (contrib == 0.0)
                    continue;

                const double prev = result.component(gp.mask);
                result.setComponent(gp.mask, prev + contrib);
            }
        }

        return result;
    }

    // Just pass null so no filter, return full product.
    inline Multivector geometricProduct(const Multivector& A, const Multivector& B) {
        return geometricProductFiltered(A, B, nullptr);
    }

}

