#pragma once

#include <bit>
#include "ga/multivector.h"
#include "ga/algebra.h"
#include "ga/signature.h"
#include "ga/basis.h"
#include "ga/ops/blade.h"

using ga::BladeMask;
using ga::Algebra;
using ga::Multivector;

namespace ga::ops {

    // Hodge dual: maps each blade to its complement blade (up to sign), using the pseudoscalar mask I_mask = (1<<dims) - 1 and geometricProductBlade
    inline Multivector dual(const Multivector& A) {
        const Algebra* alg = A.alg;
        if (!alg) {
            return A;
        }

        const int dims = alg->dimensions;
        const int bladeCount = 1 << dims;
        const BladeMask I_mask = static_cast<BladeMask>((1u << dims) - 1u);

        Multivector result(*alg);

        for (int i = 0; i < bladeCount; ++i) {
            BladeMask m = static_cast<BladeMask>(i);
            double c = A.component(m);
            if (c == 0.0)
                continue;

            // Complement mask within the n-dimensional pseudoscalar
            BladeMask comp = static_cast<BladeMask>(I_mask ^ m);

            // Blade representation with unit coefficient
            Blade a{m, +1};
            Blade b{comp, +1};

            Blade gp = geometricProductBlade(a, b, alg->signature);

            // In a well-behaved orthonormal basis, this should be ±I (or 0 in degenerate cases).
            // If gp.mask != I_mask, you may want to assert or handle specially for degenerate metrics.
            // For now, we just use the sign.
            int sign = gp.sign;

            double prev = result.component(comp);
            result.setComponent(comp, static_cast<float>(prev + c * sign));
        }

        return result;
    }

} // namespace ga::ops
