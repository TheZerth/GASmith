#pragma once
#include <bit>

namespace ga::ops {

    using Blade     = ga::Blade;
    using BladeMask = ga::BladeMask;
    using Signature = ga::Signature;

    // Robust Clifford geometric product of two basis blades.
    // Uses bitmask operations to compute:
    //   - sign from the relative ordering of basis vectors
    //   - metric contraction on overlapping axes via Signature::getSign(i)
    inline Blade geometricProductBlade(const Blade& a,
                                       const Blade& b,
                                       const Signature& sig)
    {
        // Zero handling: if either is zero, result is zero
        if (Blade::isZero(a) || Blade::isZero(b)) {
            return Blade{static_cast<BladeMask>(0), 0};
        }

        // Scalars act as identity
        if (Blade::isScalarBasis(a)) {
            return Blade{b.mask, a.sign * b.sign};
        }
        if (Blade::isScalarBasis(b)) {
            return Blade{a.mask, a.sign * b.sign};
        }

        const int dims = sig.dimensionsUsed();
        const BladeMask am = a.mask;
        const BladeMask bm = b.mask;

        // --------------------------------------------------------------------
        // 1. Sign from relative ordering (wedge part)
        // Count how many times a basis from B with smaller index would need to
        // "pass through" a basis from A.
        // This is the standard formula:
        //
        //   swaps = sum_{i in A} popcount( B & ((1 << i) - 1) )
        //
        int swapCount = 0;
        for (int i = 0; i < dims; ++i) {
            if (Blade::hasAxis(am, i)) {
                const auto lowerMask = static_cast<BladeMask>((1u << i) - 1u);
                const auto lowerInB = static_cast<BladeMask>(bm & lowerMask);
                swapCount += std::popcount(static_cast<unsigned int>(lowerInB));
            }
        }
        int sign = a.sign * b.sign * ((swapCount % 2 == 0) ? +1 : -1);

        // --------------------------------------------------------------------
        // 2. Metric contraction on overlapping axes
        //
        // Any axis present in both blades contracts:
        //   e_i e_i = g_ii
        // where g_ii = +1, -1, or 0 (null axis).
        //
        // If g_ii == 0 for any overlapping axis, the entire product is zero.
        //
        const auto overlap = static_cast<BladeMask>(am & bm);
        for (int i = 0; i < dims; ++i) {
            if (Blade::hasAxis(overlap, i)) {
                const int gii = sig.getSign(i);  // +1, -1, or 0
                if (gii == 0) {
                    // Null direction squared → zero blade
                    return Blade{static_cast<BladeMask>(0), 0};
                }
                sign *= gii;
            }
        }

        // --------------------------------------------------------------------
        // 3. Resulting basis mask
        //
        // Overlapping axes are contracted away, so the remaining basis mask is
        // simply the symmetric difference:
        //
        //   resultMask = am XOR bm
        //
        const auto resultMask = static_cast<BladeMask>(am ^ bm);

        if (sign == 0) {
            return Blade{static_cast<BladeMask>(0), 0};
        }
        if (resultMask == 0) {
            // Pure scalar result
            return Blade{static_cast<BladeMask>(0), sign};
        }

        return Blade{resultMask, sign};
    }

} // namespace ga::ops
