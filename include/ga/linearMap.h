#pragma once

#include <stdexcept>
#include <cstddef>
#include <vector>

#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/basis.h"
#include "ga/ops/wedge.h"

namespace ga {

/**
 * @brief Linear map on the vector space of an Algebra, extended to all grades
 *        by outermorphism.
 *
 * We represent L by a dense matrix m[row][col] in the chosen orthonormal basis:
 *
 *   L(e_col) = sum_row m[row][col] * e_row
 *
 * Only the first `alg->dimensions` rows/columns are used. The map is extended
 * to blades/multivectors by:
 *
 *   L(a_1 ∧ ... ∧ a_k) = L(a_1) ∧ ... ∧ L(a_k)
 *   L(α) = α
 */
struct LinearMap {
    const Algebra* alg = nullptr;
    float m[8][8];  // supports up to 8 dimensions (like the rest of GASmith)

    /// Default: no algebra, zero matrix
    LinearMap() : alg(nullptr), m{} {
        for (int r = 0; r < 8; ++r)
            for (int c = 0; c < 8; ++c)
                m[r][c] = 0.0f;
    }

    /// Construct as identity on the given Algebra
    explicit LinearMap(const Algebra& algebra) : alg(&algebra) {
        const int dims = alg->dimensions;
        for (int r = 0; r < 8; ++r) {
            for (int c = 0; c < 8; ++c) {
                m[r][c] = (r == c && r < dims) ? 1.0f : 0.0f;
            }
        }
    }

    /// Identity map factory
    static LinearMap identity(const Algebra& algebra) {
        return LinearMap(algebra);
    }

    /// Zero map factory
    static LinearMap zero(const Algebra& algebra) {
        LinearMap L;
        L.alg = &algebra;
        const int dims = algebra.dimensions;
        for (int r = 0; r < dims; ++r)
            for (int c = 0; c < dims; ++c)
                L.m[r][c] = 0.0f;
        return L;
    }

    /// Set matrix entry m[row][col] (L(e_col) component along e_row)
    void set(int row, int col, float value) {
        if (!alg)
            throw std::invalid_argument("ga::LinearMap::set: no Algebra attached");
        if (row < 0 || col < 0 || row >= alg->dimensions || col >= alg->dimensions)
            throw std::out_of_range("ga::LinearMap::set: index out of range");
        m[row][col] = value;
    }

    float get(int row, int col) const {
        if (!alg)
            throw std::invalid_argument("ga::LinearMap::get: no Algebra attached");
        if (row < 0 || col < 0 || row >= alg->dimensions || col >= alg->dimensions)
            throw std::out_of_range("ga::LinearMap::get: index out of range");
        return m[row][col];
    }

    /**
     * @brief Apply L only to the grade-1 (vector) part of v.
     *
     * Any non-vector grades are ignored. This is useful if you know v is a
     * pure vector and want a cheap application.
     */
    Multivector applyToVector(const Multivector& v) const {
        if (!alg || !v.alg || v.alg != alg) {
            throw std::invalid_argument("ga::LinearMap::applyToVector: Algebra mismatch or null");
        }

        Multivector out(*alg);
        const int dims = alg->dimensions;

        // Extract input vector components v_j
        float vcomp[8] = {0.0f};
        for (int j = 0; j < dims; ++j) {
            BladeMask mask_j = Blade::getBasis(j);
            vcomp[j] = static_cast<float>(v.component(mask_j));
        }

        // Compute w_i = sum_j m[i][j] * v_j
        for (int i = 0; i < dims; ++i) {
            float sum = 0.0f;
            for (int j = 0; j < dims; ++j) {
                sum += m[i][j] * vcomp[j];
            }
            if (sum != 0.0f) {
                out.setComponent(Blade::getBasis(i), sum);
            }
        }

        return out;
    }

    /**
     * @brief Apply the outermorphism induced by L to an arbitrary multivector A.
     *
     * Rules:
     *   - Scalars are unchanged: L(α) = α
     *   - Vectors use the linear map matrix
     *   - Higher-grade blades: L(e_{i1} ∧ ... ∧ e_{ik}) =
     *       L(e_{i1}) ∧ ... ∧ L(e_{ik})
     *   - Extended linearly to sums of blades (general multivector)
     */
    inline Multivector apply(const Multivector& A) const {
    if (!alg || !A.alg || A.alg != alg) {
        throw std::invalid_argument("ga::LinearMap::apply: Algebra mismatch or null");
    }

    using namespace ga::ops;

    const int dims = alg->dimensions;
    const std::size_t bladeCount = (1u << dims);

    // Precompute images of basis vectors: L(e_j)
    std::vector<Multivector> vecImages;
    vecImages.reserve(dims);
    for (int j = 0; j < dims; ++j) {
        Multivector ej(*alg);
        ej.setComponent(Blade::getBasis(j), 1.0f);
        vecImages.push_back(applyToVector(ej));
    }

    // Precompute images of all basis blades
    std::vector<Multivector> bladeImages;
    bladeImages.reserve(bladeCount);
    for (std::size_t i = 0; i < bladeCount; ++i) {
        bladeImages.emplace_back(*alg);
    }

    // Scalar blade (mask 0): L(1) = 1
    {
        Multivector scalarOne(*alg);
        scalarOne.setComponent(static_cast<BladeMask>(0), 1.0f);
        bladeImages[0] = scalarOne;
    }

    // For each non-scalar blade, build its image as wedge of vector images
    for (std::size_t mask = 1; mask < bladeCount; ++mask) {
        BladeMask m = static_cast<BladeMask>(mask);
        int grade = Blade::getGrade(m);
        if (grade == 1) {
            // Find which axis this mask corresponds to
            int axis = -1;
            for (int i = 0; i < dims; ++i) {
                if (m == Blade::getBasis(i)) {
                    axis = i;
                    break;
                }
            }
            if (axis < 0) {
                // Should not happen for a proper basis mask
                bladeImages[mask] = Multivector(*alg);
            } else {
                bladeImages[mask] = vecImages[axis];
            }
        } else {
            // grade >= 2: L(e_{i1} ∧ ... ∧ e_{ik}) =
            //             L(e_{i1}) ∧ ... ∧ L(e_{ik})
            BladeMask remaining = m;

            // Extract lowest set bit as first axis
            int firstAxis = -1;
            for (int i = 0; i < dims; ++i) {
                BladeMask bit = Blade::getBasis(i);
                if (Blade::hasAxis(remaining, i)) {
                    firstAxis = i;
                    remaining = static_cast<BladeMask>(remaining & ~bit);
                    break;
                }
            }

            if (firstAxis < 0) {
                bladeImages[mask] = Multivector(*alg);
                continue;
            }

            Multivector image = vecImages[firstAxis];
            if (remaining != 0) {
                Multivector restImage = bladeImages[static_cast<std::size_t>(remaining)];
                image = wedge(image, restImage);
            }

            bladeImages[mask] = image;
        }
    }

    // Now apply L by linearity:
    Multivector result(*alg);

    for (std::size_t mask = 0; mask < bladeCount; ++mask) {
        BladeMask m = static_cast<BladeMask>(mask);
        double coeff = A.component(m);
        if (coeff == 0.0)
            continue;

        const Multivector& img = bladeImages[mask];

        const int dimsLocal = alg->dimensions;
        const std::size_t N = (1u << dimsLocal);
        for (std::size_t i = 0; i < N; ++i) {
            float c = img.storage[i];
            if (c != 0.0f) {
                float prev = static_cast<float>(result.storage[i]);
                result.storage[i] = prev + static_cast<float>(coeff) * c;
            }
        }
    }

    return result;
}

};

} // namespace ga
