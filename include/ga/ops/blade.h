#pragma once
#include "../basis.h"
#include "../signature.h"

namespace ga {

    using Blade     = ga::Blade;
    using Signature = ga::Signature;

    inline Blade geometricProductBlade(const Blade& a,
                                       const Blade& b,
                                       const Signature& sig)
    {
        // Zero handling: if either is zero, result is zero
        if (Blade::isZero(a) || Blade::isZero(b)) {
            return Blade{static_cast<BladeMask>(0), 0};
        }

        // Quick scalar identity cases
        if (Blade::isScalarBasis(a)) {
            // (scalar * B) = scalar times B
            return Blade{b.mask, a.sign * b.sign};
        }
        if (Blade::isScalarBasis(b)) {
            // (A * scalar) = scalar times A
            return Blade{a.mask, a.sign * b.sign};
        }

        // Extract basis indices into arrays (sorted since your blades are canonical)
        int listA[MAX_DIMENSIONS];
        int listB[MAX_DIMENSIONS];
        int lenA = 0;
        int lenB = 0;

        for (int i = 0; i < MAX_DIMENSIONS; ++i) {
            if (Blade::hasAxis(a.mask, i)) {
                listA[lenA++] = i;
            }
            if (Blade::hasAxis(b.mask, i)) {
                listB[lenB++] = i;
            }
        }

        int sign = a.sign * b.sign;
        BladeMask resultMask = 0;

        int i = 0;
        int j = 0;

        // Core Clifford merge
        while (i < lenA && j < lenB) {
            const int ia = listA[i];
            const int jb = listB[j];

            if (ia == jb) {
                // Matching axis: e_i e_i = g_ii
                sign *= sig.getSign(ia);   // could be +1, -1, or 0
                ++i;
                ++j;
            } else if (ia < jb) {
                // Take from A directly
                resultMask |= Blade::getBasis(ia);
                ++i;
            } else {
                // jb < ia
                // Taking from B requires passing this basis through
                // the remaining ones in A -> possible sign flip
                const int remainingA = lenA - i;
                if (remainingA % 2 != 0) {
                    sign = -sign;
                }
                resultMask |= Blade::getBasis(jb);
                ++j;
            }
        }

        // Append any remaining basis from A or B
        for (; i < lenA; ++i) {
            resultMask |= Blade::getBasis(listA[i]);
        }
        for (; j < lenB; ++j) {
            resultMask |= Blade::getBasis(listB[j]);
        }

        if (sign == 0) {
            return Blade{static_cast<BladeMask>(0), 0};
        }
        if (resultMask == 0) {
            // pure scalar
            return Blade{static_cast<BladeMask>(0), sign};
        }

        return Blade{resultMask, sign};
    }

} // namespace ga
