#pragma once

#include <cmath>
#include <stdexcept>

#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/ops/geometric.h"
#include "ga/ops/involutions.h"
#include "ga/policies.h"

namespace ga {

/**
 * @brief Generic versor: an invertible multivector acting by sandwich product.
 *
 * A versor V acts on a multivector X by:
 *      X' = V X V^{-1}
 *
 * In this simple implementation we assume V ~V has a non-zero scalar part
 * and use:
 *      V^{-1} = ~V / (V ~V)_scalar
 */
class Versor {
public:
    Multivector mv;  ///< underlying multivector (should be invertible)

    Versor() = default;

    explicit Versor(const Multivector& v)
        : mv(v) {}

    /// Convenience: build from an Algebra and an already-constructed multivector.
    Versor(const Algebra& algebra, const Multivector& v)
        : mv(v)
    {
        // If v.alg is null, attach algebra; otherwise assume caller ensured consistency.
        if (!mv.alg) {
            mv.alg = &algebra;
        }
    }

    /// @return the associated Algebra (may be nullptr if mv.alg is not set).
    [[nodiscard]] const Algebra* algebra() const noexcept { return mv.alg; }

    /// @return true if this versor has a valid algebra context.
    [[nodiscard]] bool isValid() const noexcept { return mv.alg != nullptr; }

    /**
     * @brief Compute the inverse versor V^{-1}.
     *
     * Uses:
     *      V^{-1} = ~V / (V ~V)_scalar
     *
     * If no algebra is attached, or the scalar norm is (near) zero,
     * behavior is undefined (for now we just divide and let the user handle it).
     */
    [[nodiscard]] Multivector inverse() const;

    /**
     * @brief Apply the versor to a multivector X:
     *        X' = V X V^{-1}.
     */
    [[nodiscard]] Multivector apply(const Multivector& X) const;
};

// ---------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------

inline Multivector Versor::inverse() const {
    if (!mv.alg) {
        throw std::invalid_argument("ga::Versor::inverse: versor has no Algebra (mv.alg is null)");
    }

    using namespace ga::ops;

    // Reverse of the versor
    const Multivector vrev = reverse(mv);

    // For a proper versor, V * ~V is (up to metric sign) a scalar.
    const Multivector norm2_mv = geometricProduct(mv, vrev);

    // Extract scalar part (mask 0)
    const auto eps = ga::Policies::epsilon();
    if (std::fabs(s) <= eps) {
        throw std::runtime_error("ga::Versor::inverse: scalar norm (V ~V) is too close to zero");
    }

    float inv_s = 1.0f / s;

    // Scale the reversed versor
    Multivector inv = vrev;
    const int dims = mv.alg->dimensions;
    const std::size_t N = (1u << dims);
    for (std::size_t i = 0; i < N; ++i) {
        inv.storage[i] *= inv_s;
    }

    return inv;
}

inline Multivector Versor::apply(const Multivector& X) const {
    if (!mv.alg || !X.alg || mv.alg != X.alg) {
        throw std::invalid_argument("ga::Versor::apply: versor and operand must share the same Algebra");
    }

    using namespace ga::ops;
    const Multivector invV = inverse();
    const Multivector tmp  = geometricProduct(mv, X);
    return geometricProduct(tmp, invV);
}

} // namespace ga
