// C++
#pragma once
#include <cstdint>
#include <cassert>
#include <algorithm> // For std::fill
#include <cstring>   // For memset (optional, faster zeroing)

namespace ga {

struct DenseStorage {
    // Hard constraint: Max Dims = 8, so Max Elements = 2^8 = 256
    static constexpr size_t MAX_ELEMENTS = 256;

    // OPTIMIZATION 1: Use float.
    // OPTIMIZATION 2: Use a raw array (or std::array) instead of std::vector.
    // This allocates 1KB of contiguous memory directly inside the struct. ?32 bytes not 1kb
    // No heap allocation (malloc/new) ever occurs.
    float coefficients[MAX_ELEMENTS];

    // OPTIMIZATION 3: Use uint8_t for dimensions (max 255)
    uint8_t dimensions;

    explicit DenseStorage(const uint8_t dims) : dimensions(dims) {
        assert(dims <= 8 && "DenseStorage: dims too large for fixed storage");

        // Fast zeroing of memory.
        // Since we have a fixed size, we just wipe the whole 1KB.
        // It's often faster to wipe 1KB linearly than to calculate exactly
        // how much to wipe for small N.

        // OPTIMIZATION 4: Only zero the memory we'll actually use - Zur
        const size_t used_size = size_t(1) << dims;
        std::memset(coefficients, 0, used_size * sizeof(float));
    }

    // Non-const access
    float& operator[](const size_t mask) {
        // Using assert removes the check in Release builds for speed
        assert(mask < (1ULL << dimensions) && "mask out of range");
        return coefficients[mask];
    }

    // Const access
    const float& operator[](const size_t mask) const {
        assert(mask < (1ULL << dimensions) && "mask out of range");
        return coefficients[mask];
    }

    // Getter for current actual size used
    [[nodiscard]] size_t size() const {
        return size_t(1) << dimensions;
    }

    // Getter for capacity (always 256)
    [[nodiscard]] static constexpr size_t capacity() {
        return MAX_ELEMENTS;
    }
};

}