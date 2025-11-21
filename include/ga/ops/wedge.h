#pragma once

#include "ga/multivector.h"
#include "ga/ops/blade.h"
#include <bit>      // std::popcount (C++20)
#include <stdexcept>

namespace ga::ops {

inline bool keepWedgeGrade(const int gradeA, const int gradeB, const int gradeR) {
    return gradeR == gradeA + gradeB;
}

// Outer product of two multivectors
inline Multivector wedge(const Multivector& A, const Multivector& B) {
    const Multivector result{ga::ops::geometricProductFiltered(A, B, &keepWedgeGrade)};
    return result;
}

} // namespace ga::ops
