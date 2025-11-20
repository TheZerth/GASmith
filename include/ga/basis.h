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
// etc...
// A positive blade is represented by a positive coefficient. (e1^e2)
// A negative blade is represented by a negative coefficient. (-e1^e2=e2^e1)
// A zero blade is represented by a zero coefficient. (e1^e1=0)
#pragma once

#include <cstdint>
#include <array>

#define MAX_DIMENSIONS 8

namespace ga {

// You must change the BladeMask type to match your architecture. See signature.h and ensure it matches MAX_DIMENSIONS.
using BladeMask = std::uint8_t;

// Canonical Basis Blade
struct Blade {
    BladeMask mask{};
    int sign{}; // 0 for zero blade, + or - for oriented blade

// --- HELPER FUNCTIONS ---

[[nodiscard]] constexpr int getGrade(const BladeMask mask) {
    return __builtin_popcount(mask);
}

[[nodiscard]] constexpr bool hasAxis(BladeMask mask, int i) {
    bool result = false;
    if (i >=0 && i < 8) {
        result = (mask & (BladeMask(1u) << i)) != 0; // does mask contain bit i?
    }
    return result;
}

[[nodiscard]] constexpr Blade getBasis(int axisIndex) {
    return (BladeMask(1u) << axisIndex, 1);
}

[[nodiscard]] constexpr int highestAxis(BladeMask mask) {
    return 8 - __builtin_clz(mask); // clz = count leading zeros
}

[[nodiscard]] constexpr bool doesOverlap(BladeMask a, BladeMask b) {
    return (a & b) != 0;
}

[[nodiscard]] constexpr BladeMask addAxis(BladeMask mask, int axisIndex) {
    return mask | (BladeMask(1u) << axisIndex);
}

[[nodiscard]] constexpr BladeMask removeAxis(BladeMask mask, int axisIndex) {
    return mask & ~(BladeMask(1u) << axisIndex);
}

[[nodiscard]] constexpr BladeMask toggleAxis(BladeMask mask, int axisIndex) {
    return mask ^ (BladeMask(1u) << axisIndex);
}

// --- CONSTRUCTORS ---

constexpr Blade() = default;

constexpr Blade(BladeMask mask, int sign) : mask(mask), sign(sign) {};

constexpr Blade makeBlade(const int* basis, int numBasis) {
    if (numBasis <= 0) {
        return Blade(BladeMask(0), 1); // Scalar
    }
    if (numBasis > MAX_DIMENSIONS) {
        return Blade(BladeMask(0), 0); // Treat as zero blade for safety
    }
    if (!basis) {
        return Blade(BladeMask(0), 0);
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

    // selection sort + swap count
    for (int i = 0; i + 1 < numBasis; ++i) {
        int minIdx = i;
        for (int j = i + 1; j < numBasis; ++j) {
            if (tempBasis[static_cast<std::size_t>(j)]
                < tempBasis[static_cast<std::size_t>(minIdx)]) {
                minIdx = j;
                }
        }
        if (minIdx != i) {
            std::swap(tempBasis[static_cast<std::size_t>(i)],
                      tempBasis[static_cast<std::size_t>(minIdx)]);
            ++swaps;
        }
    }

    // parity of swaps → sign
    result.sign = (swaps % 2 == 0) ? +1 : -1;
    uint8_t mask = 0;
    for (int i = 0; i < numBasis; ++i) {
        const int idx = basis[static_cast<std::size_t>(i)];
        // optional: bounds check 0 <= idx < MAX_DIMENSIONS
        mask |= BladeMask(1u << idx);
    }
    result.mask = mask;
    return result;
}

};

}