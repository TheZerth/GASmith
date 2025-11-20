// --- SIMPLE ---
// A signature defines the rules of measurement of our space.
// This tells our algebra how distances and angles behave.
// Signatures are composed of three components: Positive(P), Negative(N), and Zero(R).
// Each axis is assigned to one of these components. This defines how the axis behaves.
// A P axis, acts like a standard number line. A length squared is positive.
// An Q axis, acts like an imaginary number line. A length squared is negative.
// An R axis, acts like a null line. A length squared is zero.
// When defining a signature, we specify the number of P, N, and R axes.
// For example, a signature of (3, 1, 0) defines a space with 3 P axes, 1 N axis, and 0 R axes.
//
// 3d Euclidean Space is represented by the signature (3, 0, 0). 3 axis (3D)
// STA uses the signature (1, 3, 0). 4 axis (1 time + 3 space)
// PGA uses the signature (3, 0, 1). 4 axis (3D space + 1 point at infinity)
// CGA uses the signature (4, 1, 0). 5 axis
// Etc... many more signatures are possible depending on the space being described.

// --- COMPLEX ---
// A signature is represented as a triplet (P, N, R).
// This is a metric signature specifying the linear mapping of g(i,j)
// ... TO DO

#pragma once
#include <array>
#include <cstdint>
#include <stdexcept>

namespace ga {

// Maximum number of dimensions supported by the algebra. (Capacity, not always fully used)
static constexpr int MAX_DIMENSIONS = 8; // Almost all use cases are in 3,4 or 5 dimensions. Limiting to 8 for extra support.
// If performing clifford CL(p,q) analysis, exotic spaces analysis or quantum information analysis, adjust MAX_DIMENSIONS as needed.
// Increasing will have significant impact on performance the signature is used to generate the blades and the algebras.
// Blades are defined by 2^n where n is # of dimensions.
// For 3D Euclidean: 2^3, 8 coefficients, easy
// For 3D CGA: 2^5, 32, doable
// For Conformal Space Time: 2^6, 64, medium
// Exotic Space: 2^16, 65536 coefficients, HARD
// DO NOT EXCEED 16, our dense storage system will not work beyond 16D. Need to change to sparse storage if above 16D.

using Metric = std::array<int, MAX_DIMENSIONS>;
using Mask = std::array<bool, MAX_DIMENSIONS>;

struct Signature {
private:
    // --- VARIABLES ---
    // Data
    int p_{}; // # of positive axis
    int q_{}; // # of negative axis
    int r_{}; // # of null axis
    Metric metric_{}; // g_ij metric coefficient, assuming right angle basis vectors so single dimension vector. TO DO expand to non-diagonal metric.
    bool isRightHanded_{true}; // Orientation flag; A convention must be agreed upon which orientation is positive.

    // Helpers
    int dimensionsUsed_{};

public:
    // --- FUNCTIONS ---
    // Getters
    [[nodiscard]] constexpr int p() const { return p_; }
    [[nodiscard]] constexpr int q() const { return q_; }
    [[nodiscard]] constexpr int r() const { return r_; }
    [[nodiscard]] constexpr Metric metric() const { return metric_; }

    // Metric Lookup
    [[nodiscard]] constexpr int metricLookup(const int i, const int j) const {
        if (i < MAX_DIMENSIONS && j < MAX_DIMENSIONS && i >= 0 && j >= 0) {
            return (i==j) ? metric_[i] : 0;
        }
        return -2;
    }  // return the value of g at index i and j. Assumes an orthogonal basis. TO DO expand to non-diagonal metric

    // Helpers
    [[nodiscard]] constexpr int dimensionsUsed() const { return dimensionsUsed_; }
    [[nodiscard]] constexpr int getSign(const int i = 0) const { return metricLookup(i,i); }
    [[nodiscard]] constexpr bool isRightHanded() const { return isRightHanded_; }
    [[nodiscard]] constexpr bool isLeftHanded() const { return !isRightHanded_; }
    [[nodiscard]] constexpr bool isPos(const int i) const { return (metricLookup(i,i) == 1); }
    [[nodiscard]] constexpr bool isNeg(const int i) const { return (metricLookup(i,i) == -1); }
    [[nodiscard]] constexpr bool isZero(const int i) const { return (metricLookup(i,i) == 0); }
    [[nodiscard]] constexpr bool isDegenerate() const { return (r_ > 0); } // the signature is degenerate if the algebra contains a null axis (math jargon)

    // Default constructor, no dimensions.
    constexpr Signature() = default;

    // Construct a signature given known dimensions.
    constexpr Signature(const int p, const int q, const int r, const bool isRightHanded)
                        : p_{p}, q_{q}, r_{r}, isRightHanded_{isRightHanded} {
        if (!buildMetric(p, q, r)) {
            throw std::invalid_argument("Failed to construct signature. buildMetric() failed, exceeds max dimensions.");
        }
    }

    // Construct a signature given known p, q and r masks.
    constexpr Signature(const Mask pMask, const Mask qMask, const Mask rMask,
                        const bool isRightHanded) : isRightHanded_{isRightHanded} {
        if (!buildMetric(pMask, qMask, rMask)) {
            throw std::invalid_argument("Failed to construct signature. buildMetric() failed, overlapping mask.");
        }
    }

    // Construct a signature given a diagonal metric. TO DO expand to non-diagnonal metric
    constexpr Signature(const Metric &metric, const int axisCount, const bool isRightHanded) : metric_{metric}, isRightHanded_{isRightHanded} {
        if (axisCount > MAX_DIMENSIONS || axisCount < 0) {
            throw std::invalid_argument("Failed to construct signature. axisCount out of range.");
        }
        extractMetric(metric, axisCount);
    }

    // Takes a metric and extracts values to p q and r; TO DO extend to non-diagonal metric
    constexpr void extractMetric(const Metric &metric, const int axisCount) {
        // Build temp axis counts
        int p = 0;
        int q = 0;
        int r = 0;
        // Count axis
        for (int i = 0; i < axisCount; i++) {
            switch (metric[i]) {
                case 1: ++p; break;
                case -1: ++q; break;
                default: ++r; break;
            }
        }
        // Assign extracted values
        p_ = p; q_ = q; r_ = r;
        dimensionsUsed_ = p_ + q_ + r_;
    }

    // Takes p q and r and constructs a diagonal metric; TO DO extend to non-diagonal metric
    constexpr bool buildMetric(const int p, const int q, const int r) {
        // Ensure does not exceed max
        if (p + q + r > MAX_DIMENSIONS)
            return false;
        // Create blank metric and iterator
        Metric metric = {};
        int i = 0;
        // Assign positive axis
        for (int n = 0; n < p; ++n) {
            metric[i++] = 1;
        }
        // Assign negative axis
        for (int n = 0; n < q; ++n) {
            metric[i++] = -1;
        }
        // Assign null axis
        for (int n = 0; n < r; ++n) {
            metric[i++] = 0;
        }
        // Assign metric
        metric_ = metric;
        // Assign extracted values
        p_ = p; q_ = q; r_ = r;
        dimensionsUsed_ = p_ + q_ + r_;
        return true;
    }

    // Takes p q r masks and constructs a diagonal metric; TO DO extend to non-diagonal metric
    constexpr bool buildMetric(const Mask pMask, const Mask qMask, const Mask rMask) {
        // Ensure no mask overlap
        if (!validateMasks(pMask, qMask, rMask))
            return false;
        // Create blank metric
        Metric metric = {};
        int n = 0;

        dimensionsUsed_ = 0;
        p_ = 0, q_ = 0, r_ = 0;
        for (int i = 0; i < MAX_DIMENSIONS; ++i) {

            // Assign positive axis
            if (pMask[i]) {
                metric[i] = 1;
                ++n, ++p_, ++dimensionsUsed_;
            }

            // Assign negative axis
            if (qMask[i]) {
                metric[i] = -1;
                ++n, ++p_, ++dimensionsUsed_;
            }

            // Assign null axis
            if (rMask[i]) {
                metric[i] = 0;
                ++n, ++p_, ++dimensionsUsed_;
            }
        }

        // Assign metric
        metric_ = metric;
        return true;
    }

    // Validates that p q and r masks do not overlap
    static constexpr bool validateMasks(const Mask pMask, const Mask qMask, const Mask rMask) {
        for (int i = 0; i < MAX_DIMENSIONS; ++i) {
            if ((pMask[i] && qMask[i]) || (pMask[i] && rMask[i]) || (qMask[i] && rMask[i])) {
                return false;
            }
        }
        return true;
    }
};

} // namespace ga