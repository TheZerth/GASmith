#pragma once

#include <stdexcept>

#include "blade.h"
#include "ga/algebra.h"
#include "ga/multivector.h"

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
                                           const GradeFilterFn keep) {
        const Algebra *alg = A.alg;

        if (!alg || !B.alg || alg != B.alg) {
            throw std::invalid_argument(
                    "ga::ops::geometricProductFiltered: Multivectors must share the same Algebra");
        }

        const int dims = alg->dimensions;
        const int bladeCount = 1 << dims;

        Multivector result(*alg);

        // Pre-compute grades if filtering is enabled
        const bool useFilter = (keep != nullptr);

        for (int i = 0; i < bladeCount; ++i) {
            const auto maskA = static_cast<BladeMask>(i);
            const double coeffA = A.component(maskA);
            if (coeffA == 0.0)
                continue;

            const int gradeA = useFilter ? ga::Blade::getGrade(maskA) : 0;

            for (int j = 0; j < bladeCount; ++j) {
                const auto maskB = static_cast<BladeMask>(j);
                const double coeffB = B.component(maskB);
                if (coeffB == 0.0)
                    continue;

                const int gradeB = useFilter ? ga::Blade::getGrade(maskB) : 0;

                const Blade gp = geometricProductBlade(
                        Blade{maskA, +1},
                        Blade{maskB, +1},
                        alg->signature
                );

                if (Blade::isZero(gp))
                    continue;

                if (useFilter) {
                    const int gradeR = ga::Blade::getGrade(gp.mask);
                    if (!keep(gradeA, gradeB, gradeR))
                        continue;
                }

                const double contrib = coeffA * coeffB * static_cast<double>(gp.sign);
                result.setComponent(gp.mask, result.component(gp.mask) + contrib);
            }
        }

        return result;
    }


        // Just pass null so no filter, return full product.
        inline Multivector geometricProduct(const Multivector& A, const Multivector& B) {
            return geometricProductFiltered(A, B, nullptr);
        }

}


