#pragma once

#include "ga/multivector.h"
#include "ga/ops/blade.h"
#include <bit>      // std::popcount (C++20)
#include <stdexcept>

namespace ga::ops {

// Helper: grade of a basis blade given its mask
inline int gradeOfMask(BladeMask m) {
    return std::popcount(static_cast<unsigned int>(m));
}

// Outer product of two multivectors
inline Multivector wedge(const Multivector& A, const Multivector& B) {
    const Algebra* algA = A.alg;
    const Algebra* algB = B.alg;

    if (!algA || !algB || algA != algB) {
        throw std::invalid_argument("ga::ops::wedge: Multivectors must share the same Algebra");
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

        int gradeA = gradeOfMask(maskA);

        for (int j = 0; j < bladeCount; ++j) {
            BladeMask maskB = static_cast<BladeMask>(j);
            double coeffB = B.component(maskB);
            if (coeffB == 0.0)
                continue;

            int gradeB = gradeOfMask(maskB);

            // Build unit coefficient blades for geometricProductBlade
            Blade a{maskA, +1};
            Blade b{maskB, +1};

            Blade gp = geometricProductBlade(a, b, alg->signature);

            int resultGrade = gradeOfMask(gp.mask);

            // Keep only the "outer" grade
            if (resultGrade != gradeA + gradeB)
                continue;

            double contrib = coeffA * coeffB * static_cast<double>(gp.sign);
            if (contrib == 0.0)
                continue;

            double prev = result.component(gp.mask);
            result.setComponent(gp.mask, prev + contrib);
        }
    }

    return result;
}

} // namespace ga::ops
