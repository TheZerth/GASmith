#pragma once

#include <cmath>
#include <stdexcept>

#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/versor.h"
#include "ga/ops/geometric.h"
#include "ga/ops/involutions.h"
#include "ga/ops/wedge.h"
#include "ga/ops/inner.h"
#include "ga/policies.h"

namespace ga {

/**
 * @brief Rotor: an even versor representing a (metric-dependent) rotation / Lorentz transformation.
 *
 * In Euclidean spaces, a rotor R acts on a vector (or multivector) X by:
 *      X' = R X ~R
 *
 * and satisfies R ~R = 1 (unit rotor).
 *
 * This implementation does not enforce "even grade" structurally, but the user
 * is expected to construct rotors via fromPlaneAngle or other even-grade constructions.
 */
class Rotor {
public:
    Multivector mv;  ///< underlying multivector, expected to be a unit versor

    Rotor() = default;

    explicit Rotor(const Multivector& R)
        : mv(R) {}

    /// @return the associated Algebra (may be nullptr if mv.alg is not set).
    const Algebra* algebra() const noexcept { return mv.alg; }

    /// @return true if rotor has a valid algebra.
    bool isValid() const noexcept { return mv.alg != nullptr; }

    /// @return underlying multivector (const)
    const Multivector& value() const noexcept { return mv; }

    /// @return underlying multivector (mutable)
    Multivector& value() noexcept { return mv; }

    /**
     * @brief Normalize the rotor so that R ~R = 1 (up to numerical precision).
     *
     * Uses:
     *      norm2 = R ~R (scalar)
     *      R <- R / sqrt(norm2)
     */
    void normalize();

    /**
     * @brief Apply the rotor to a multivector X:
     *        X' = R X ~R.
     *
     * In Euclidean 3D, this is a proper rotation. In other signatures it becomes
     * a metric-appropriate Lorentz-like transformation.
     */
    Multivector apply(const Multivector& X) const;

    /**
     * @brief Construct a rotor from a plane (bivector) and angle.
     *
     * In Euclidean 3D, given a unit bivector B describing a plane of rotation,
     * the rotor is:
     *      R = cos(theta/2) - B sin(theta/2)
     *
     * This function assumes B is (approximately) unit with respect to the metric.
     * You can normalize B beforehand if needed.
     */
    static Rotor fromBivectorAngle(const Multivector& B, float theta);

    /**
     * @brief Construct a rotor from two vectors a, b defining the rotation plane
     *        and angle.
     *
     * We form the bivector B = a ∧ b and then construct:
     *      R = cos(theta/2) - \hat{B} sin(theta/2)
     *
     * where \hat{B} is B normalized. For now this routine does a simple
     * magnitude estimate; for production use you may wish to provide a
     * more robust metric-aware normalization.
     */
    static Rotor fromPlaneAngle(const Multivector& a,
                                const Multivector& b,
                                float theta);
};

// ---------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------

inline void Rotor::normalize() {
    if (!mv.alg) {
        throw std::invalid_argument("ga::Rotor::normalize: rotor has no Algebra (mv.alg is null)");
    }

    using namespace ga::ops;

    Multivector rrev = reverse(mv);
    Multivector norm2_mv = geometricProduct(mv, rrev);

    float s = static_cast<float>(norm2_mv.component(static_cast<ga::BladeMask>(0)));
    const auto eps = ga::Policies::epsilon();
    if (std::fabs(s) <= eps) {
        throw std::runtime_error("ga::Rotor::normalize: rotor norm^2 is too close to zero");
    }

    float inv_sqrt = 1.0f / std::sqrt(std::fabs(s));

    const int dims = mv.alg->dimensions;
    const std::size_t N = (1u << dims);
    for (std::size_t i = 0; i < N; ++i) {
        mv.storage[i] *= inv_sqrt;
    }
}

inline Multivector Rotor::apply(const Multivector& X) const {
    if (!mv.alg || !X.alg || mv.alg != X.alg) {
        throw std::invalid_argument("ga::Rotor::apply: rotor and operand must share the same Algebra");
    }

    using namespace ga::ops;

    Multivector rrev = reverse(mv);
    Multivector tmp  = geometricProduct(mv, X);
    return geometricProduct(tmp, rrev);
}

inline Rotor Rotor::fromBivectorAngle(const Multivector& B, float theta) {
    if (!B.alg) {
        throw std::invalid_argument("ga::Rotor::fromBivectorAngle: bivector has no Algebra");
    }

    const Algebra* alg = B.alg;
    const int dims = alg->dimensions;
    const std::size_t N = (1u << dims);

    Multivector R(*alg);

    // Scalar part: cos(theta/2)
    float c = std::cos(theta * 0.5f);
    R.setComponent(static_cast<ga::BladeMask>(0), c);

    // Bivector part: -sin(theta/2) * B (assumes B is unit)
    float s = std::sin(theta * 0.5f);
    for (std::size_t i = 0; i < N; ++i) {
        float bcoef = B.storage[i];
        if (bcoef != 0.0f) {
            R.storage[i] += -s * bcoef;
        }
    }

    Rotor rot(R);
    rot.normalize();
    return rot;
}

inline Rotor Rotor::fromPlaneAngle(const Multivector& a,
                                   const Multivector& b,
                                   float theta)
{
    if (!a.alg || !b.alg || a.alg != b.alg) {
        throw std::invalid_argument("ga::Rotor::fromPlaneAngle: a and b must share the same Algebra");
    }

    using namespace ga::ops;

    const Algebra* alg = a.alg;

    // Bivector B = a ∧ b
    Multivector B = wedge(a, b);

    // Metric-aware magnitude of B using the inner product: for a pure bivector, B ⋅ B is scalar.
    Multivector BB = inner(B, B);
    float norm2 = static_cast<float>(BB.component(static_cast<ga::BladeMask>(0)));
    const auto eps = ga::Policies::epsilon();
    if (std::fabs(norm2) <= eps) {
        throw std::runtime_error(
            "ga::Rotor::fromPlaneAngle: a ∧ b has zero (or near-zero) norm (no well-defined plane)");
    }

    float inv_mag = 1.0f / std::sqrt(std::fabs(norm2));

    // Normalize B using the metric-aware magnitude
    const int dims = alg->dimensions;
    const std::size_t N = (1u << dims);

    // Normalize B (combinatorial, not fully metric-aware, but good enough for now)
    for (std::size_t i = 0; i < N; ++i) {
        B.storage[i] *= inv_mag;
    }

    return fromBivectorAngle(B, theta);
}

} // namespace ga
