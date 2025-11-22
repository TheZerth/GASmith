#pragma once

#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/ops/wedge.h"

namespace ga {

    struct LinearMap {
        const Algebra* alg;
        float m[8][8]; // row-major, only dims used

        explicit LinearMap(const Algebra& algebra);

        // L(e_j) = sum_i m[i][j] e_i
        Multivector applyToVector(const Multivector& v) const;

        // O(2^n) outermorphism: extend to whole multivector
        Multivector apply(const Multivector& A) const;
    };

} // namespace ga
