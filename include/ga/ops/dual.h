#pragma once
#include "ga/ops/geometric.h"
#include "ga/multivector.h"
#include "ga/algebra.h"

namespace ga::ops {

    using ga::BladeMask;
    using ga::Algebra;
    using ga::Multivector;
    using ga::Blade;

    // Hodge dual: maps each blade to its complement blade (up to sign), using the pseudoscalar mask I_mask = (1<<dims) - 1 and geometricProductBlade
    inline Multivector dual(const Multivector& A) {
        const Algebra* alg = A.alg;
        if (!alg) {
            return A;
        }

        const int dims = alg->dimensions;
        const int bladeCount = 1 << dims;
        const auto I_mask = static_cast<BladeMask>((1u << dims) - 1u);

        Multivector result(*alg);

        for (int i = 0; i < bladeCount; ++i) {
            const auto m = static_cast<BladeMask>(i);
            const double c = A.component(m);
            if (c == 0.0)
                continue;

            // Complement mask within the n-dimensional pseudoscalar
            const auto comp = static_cast<BladeMask>(I_mask ^ m);

            // Make degeneracy handling explicit: if the product does not yield the pseudoscalar
            // or has zero sign, we treat the dual as undefined for this blade and skip it.
            const auto gp = ga::ops::geometricProductBlade(
                    Blade{m, +1},
                    Blade{comp, +1},
                    alg->signature
            );

            if (gp.sign == 0 || gp.mask != I_mask) {
                // Degenerate or ill-defined dual for this blade; skip contribution.
                continue;
            }

            const int sign = gp.sign;

            result.setComponent(comp, result.component(comp) + c * sign);
        }

        return result;
    }

} // namespace ga::ops
