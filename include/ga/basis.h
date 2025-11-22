// --- SIMPLE ---
// A basis defines the building blocks of the algebra
// Typically represented e1, e2, e3... eN
// These basis are our axis. Ex: x = e1, y = e2, z = e3
//
// We can combine basis to construct blades
// Blades are constructed by taking the outer product of one or more basis
// e1 (1-vector) - a vector along axis e1
// e2^e3 (2-vector) - an oriented plane along axis e2 e3
// e1^e2^e3 (3-vector) - an oriented cube along axis e1 e2 e3
// etc...
//
// Blades are represented using bitmasks and a coefficent for the orientation.
// In a bitmask every bit represents a basis. And 0 or 1 encode if the basis is present or not.
// Ex: 00000001 = e1, 00000010 = e2, 00000100 = e3
// 00001010 = e2^e4, 00010100 = e3^e5
// 00101010 = e1^e2^e3
// etc... Note that the index is -1 so basis index 0 is e1.
// A positive blade is represented by a positive coefficient. (e1^e2)
// A negative blade is represented by a negative coefficient. (-e1^e2=e2^e1)
// A zero blade is represented by a zero coefficient. (e1^e1=0)
//
// Usage Note
// A blade of {static_cast<BladeMask>(0), 1}; represents the unit scalar basis or "1"
// A blade of {static_cast<BladeMask>(0), 0}; represents the zero blade or a wedge collapse. Equivalent to 0.
#pragma once
#include <cstdint>
#include <array>

static constexpr int MAX_DIMENSIONS = 8;

namespace ga {

    // You must change the BladeMask type to match your architecture. See signature.h and ensure it matches MAX_DIMENSIONS.
    using BladeMask = std::uint8_t;

    // Canonical Basis Blade
    struct Blade {
        BladeMask mask{};
        int sign{}; // 0 for zero blade, + or - for oriented blade


        // --- HELPER FUNCTIONS ---
        [[nodiscard]] static constexpr int getGrade(const BladeMask mask) {
            // BladeMask is uint8_t-compatible; popcount wants a wider integer
            return std::popcount(static_cast<unsigned int>(mask));
        }

        [[nodiscard]] static constexpr bool hasAxis(const BladeMask mask, const int i) {
            bool result = false;
            if (i >=0 && i < 8) {
                result = (mask & (static_cast<BladeMask>(1u) << i)) != 0; // does mask contain bit i?
            }
            return result;
        }

        // Return a basis vector index as a bitmask
        [[nodiscard]] static constexpr BladeMask getBasis(const int axisIndex) {
            if (axisIndex < 0 || axisIndex >= MAX_DIMENSIONS) return 0;
            return static_cast<BladeMask>(1u) << axisIndex; // Shift bit over by index and return
        }

        [[nodiscard]] static constexpr int highestAxis(const BladeMask mask) {
            if (mask == 0) return -1;

            int highest = -1;
            for (int i = MAX_DIMENSIONS - 1; i >= 0; --i) {
                if (mask & (static_cast<BladeMask>(1u) << i)) {
                    highest = i;
                    break;
                }
            }
            return highest;
        }

        [[nodiscard]] static constexpr bool doesOverlap(const BladeMask a, const BladeMask b) {
            return (a & b) != 0;
        }

        [[nodiscard]] static constexpr bool isZero(const Blade blade) {
            return blade.sign == 0;
        }

        [[nodiscard]] static constexpr bool isScalarBasis(const Blade blade) {
            return (blade.mask == 0) && (blade.sign != 0);
        }

        [[nodiscard]] static constexpr BladeMask addAxis(const BladeMask mask, const int axisIndex) {
            return mask | (static_cast<BladeMask>(1u) << axisIndex);
        }

        [[nodiscard]] static constexpr BladeMask removeAxis(const BladeMask mask, const int axisIndex) {
            return mask & ~(static_cast<BladeMask>(1u) << axisIndex);
        }

        [[nodiscard]] static constexpr BladeMask toggleAxis(const BladeMask mask, const int axisIndex) {
            return mask ^ (static_cast<BladeMask>(1u) << axisIndex);
        }

    // --- CONSTRUCTORS ---

        constexpr Blade() = default;

        constexpr Blade(const BladeMask mask, const int sign) : mask(mask), sign(sign) {};

        static constexpr Blade makeBlade(const int* basis, const int numBasis) {
            if (numBasis <= 0) {
                return {static_cast<BladeMask>(0), 1}; // Scalar
            }
            if (numBasis > MAX_DIMENSIONS) {
                return {static_cast<BladeMask>(0), 0}; // Treat as zero blade for safety
            }
            if (!basis) {
                return {static_cast<BladeMask>(0), 0};
            }

            // Create empty blade
            Blade result{};
            // Create buffer for sorting
            std::array<int, MAX_DIMENSIONS> tempBasis{};

            // copy into local buffer
            for (int i = 0; i < numBasis; ++i) {
                tempBasis[i] = basis[i];
            }

            int swaps = 0;

            // bubble sort + swap count
            for (int i = 0; i < numBasis - 1; ++i) {
                for (int j = 0; j < numBasis - 1 - i; ++j) {
                    if (tempBasis[j] > tempBasis[j + 1]) {
                        std::swap(tempBasis[j], tempBasis[j + 1]);
                        ++swaps;
                    } else if (tempBasis[j] == tempBasis[j + 1]) {
                        // Duplicate axis → wedge is zero
                        return Blade{static_cast<BladeMask>(0), 0};
                    }
                }
            }

            // parity of swaps → sign
            result.sign = (swaps % 2 == 0) ? +1 : -1;
            uint8_t mask = 0;
            for (int i = 0; i < numBasis; ++i) {
                const int idx = tempBasis[i];
                // optional: bounds check 0 <= idx < MAX_DIMENSIONS
                mask |= static_cast<BladeMask>(1u << idx);
            }
            result.mask = mask;
            return result;
        }

        static constexpr Blade combineBlade(const Blade a, const Blade b) {
            // Zero blade in, zero out
            if (isZero(a) || isZero(b)) {
                return Blade{static_cast<BladeMask>(0), 0};
            }

            // Scalar basis identity: 1 ^ B = B,  B ^ 1 = B
            if (isScalarBasis(a)) {
                return Blade{b.mask, a.sign * b.sign};
            }
            if (isScalarBasis(b)) {
                return Blade{a.mask, a.sign * b.sign};
            }

            // Overlap -> wedge = 0
            if (doesOverlap(a.mask, b.mask)) {
                return Blade{static_cast<BladeMask>(0), 0};
            }

            const BladeMask resultMask = a.mask ^ b.mask; // Blades are already sorted so if no overlaps just superimpose.

            // Compute parity: how many "swaps" needed to move b's axes past a's axes.
            // This is the same as your wedgeParity(..) function if you have one.
            int swaps = 0;
            for (int j = 0; j < MAX_DIMENSIONS; ++j) {
                if (!hasAxis(b.mask, j)) continue;
                for (int i = j + 1; i < MAX_DIMENSIONS; ++i) {
                    if (hasAxis(a.mask, i)) {
                        ++swaps;
                    }
                }
            }

            const int paritySign = (swaps % 2 == 0) ? +1 : -1;
            const int resultSign = a.sign * b.sign * paritySign;

            return Blade{resultMask, resultSign};
        }

    };

}