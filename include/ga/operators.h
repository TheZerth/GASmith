#pragma once

#include "ga/multivector.h"
#include "ga/ops/geometric.h"
#include "ga/ops/wedge.h"
#include "ga/ops/inner.h"

namespace ga {

// Geometric product: A * B
inline Multivector operator*(const Multivector& A, const Multivector& B) {
    return ga::ops::geometricProduct(A, B);
}

// Wedge (outer) product: A ^ B
inline Multivector operator^(const Multivector& A, const Multivector& B) {
    return ga::ops::wedge(A, B);
}

// Hestenes inner product: A & B
inline Multivector operator&(const Multivector& A, const Multivector& B) {
    return ga::ops::inner(A, B);
}

// Left/right contraction: A << B and A >> B
inline Multivector operator<<(const Multivector& A, const Multivector& B) {
    return ga::ops::leftContraction(A, B);
}

inline Multivector operator>>(const Multivector& A, const Multivector& B) {
    return ga::ops::rightContraction(A, B);
}

inline Multivector operator~(const Multivector& A) {
    return ga::ops::reverse(A);
}

inline Multivector operator!(const Multivector& A) {
    return ga::ops::cliffordConjugate(A);
}

inline Multivector dual(const Multivector& A) {
    return ga::ops::dual(A);
}

inline std::ostream& operator<<(std::ostream& os, const Multivector& A) {
    if (!A.alg) {
        return os << "Multivector{<null algebra>}";
    }
    const int dims = A.alg->dimensions;
    const std::size_t N = (1u << dims);

    bool first = true;
    for (std::size_t i = 0; i < N; ++i) {
        float c = A.storage[i];
        if (c == 0.0f) continue;

        if (!first) os << " + ";
        first = false;

        os << c;

        if (i != 0) {
            os << "e";
            BladeMask m = static_cast<BladeMask>(i);
            bool firstAxis = true;
            for (int axis = 0; axis < dims; ++axis) {
                if (Blade::hasAxis(m, axis)) {
                    if (!firstAxis) os << "";
                    os << (axis + 1);
                    firstAxis = false;
                }
            }
        }
    }

    if (first) {
        os << "0";
    }

    return os;
}

} // namespace ga
