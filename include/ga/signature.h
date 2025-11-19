//
// Created by zerth on 10/17/25.
//

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

static constexpr int ALGEBRA_DIMENSIONS = 4;

using Metric = std::array<int, ALGEBRA_DIMENSIONS>;
using Mask = std::array<bool, ALGEBRA_DIMENSIONS>;

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
    [[nodiscard]] constexpr int metricLookup(const int i, const int j) const { return (i==j) ? metric_[i] : 0; }  // return the value of g at index i and j. Assumes an orthogonal basis. TO DO expand to non-diagonal metric

    // Helpers
    [[nodiscard]] constexpr int dimensionsUsed() const { return dimensionsUsed_; }
    [[nodiscard]] constexpr int getSign(int i = 0) const { return metricLookup(i,i); }
    [[nodiscard]] constexpr bool isRightHanded() const { return isRightHanded_; }
    [[nodiscard]] constexpr bool isLeftHanded() const { return !isRightHanded_; }
    [[nodiscard]] constexpr bool isPos(int i) const { return (metricLookup(i,i) == 1); }
    [[nodiscard]] constexpr bool isNeg(int i) const { return (metricLookup(i,i) == -1); }
    [[nodiscard]] constexpr bool isZero(int i) const { return (metricLookup(i,i) == 0); }
    [[nodiscard]] constexpr bool isDegenerate() const { return (r_ > 0); } // the signature is degenerate if the algebra contains a null axis (math jargon)

    // Default constructor, no dimensions.
    constexpr Signature() = default;

    // Construct a signature given known dimensions.
    constexpr Signature(const int p, const int q, const int r, const bool isRightHanded)
                        : p_{p}, q_{q}, r_{r}, isRightHanded_{isRightHanded} {
        if (!buildMetric(p, q, r)) {
            throw std::invalid_argument("Failed to construct signature. buildMetric() failed, p q r sum does not equal algebra dimensions.");
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
    constexpr Signature(const Metric metric, const bool isRightHanded) : metric_{metric}, isRightHanded_{isRightHanded} {
        extractMetric(metric);
    }

    // Takes a metric and extracts values to p q and r; TO DO extend to non-diagonal metric
    constexpr void extractMetric(const Metric metric) {
        // Build temp axis counts
        int p = 0;
        int q = 0;
        int r = 0;
        // Count axis
        for (int i : metric) {
            switch (i) {
                case 1: p++; break;
                case -1: q++; break;
                default: r++; break;
            }
        }
        // Assign extracted values
        p_ = p; q_ = q; r_ = r;
        dimensionsUsed_ = p_ + q_ + r_;
    }

    // Takes p q and r and constructs a diagonal metric; TO DO extend to non-diagonal metric
    constexpr bool buildMetric(const int p, const int q, const int r) {
        // Ensure no mask overlap
        if (p + q + r != ALGEBRA_DIMENSIONS)
            return false;
        // Create blank metric and iterator
        Metric metric = {};
        int i = 0;
        // Assign positive axis
        for (int n = 0; n < p; n++) {
            metric[i] = 1;
            i++;
        }
        // Assign negative axis
        for (int n = 0; n < q; n++) {
            metric[i] = -1;
            i++;
        }
        // Assign null axis
        for (int n = 0; n < r; n++) {
            metric[i] = 0;
            i++;
        }
        // Assign metric
        metric_ = metric;
        extractMetric(metric_);
        return true;
    }

    // Takes p q r masks and constructs a diagonal metric; TO DO extend to non-diagonal metric
    constexpr bool buildMetric(const Mask pMask, const Mask qMask, const Mask rMask) {
        // Ensure no mask overlap
        if (!validateMasks(pMask, qMask, rMask))
            return false;
        // Create blank metric
        Metric metric = {};
        // Assign positive axis
        for (int i = 0; i < ALGEBRA_DIMENSIONS; i++) {
            if (pMask[i]) {
                metric[i] = 1;
            }
        }
        // Assign negative axis
        for (int i = 0; i < ALGEBRA_DIMENSIONS; i++) {
            if (qMask[i]) {
                metric[i] = -1;
            }
        }
        // Assign null axis
        for (int i = 0; i < ALGEBRA_DIMENSIONS; i++) {
            if (rMask[i]) {
                metric[i] = 0;
            }
        }
        // Assign metric
        metric_ = metric;
        extractMetric(metric_);
        return true;
    }

    // Validates that p q and r masks do not overlap
    static constexpr bool validateMasks(const Mask pMask, const Mask qMask, const Mask rMask) {
        for (int i = 0; i < ALGEBRA_DIMENSIONS; i++) {
            if (pMask[i] && qMask[i] == true) {
                return false;
            }
            if (pMask[i] && rMask[i] == true) {
                return false;
            }
            if (qMask[i] && rMask[i] == true) {
                return false;
            }
        }
        return true;
    }
};

} // namespace ga