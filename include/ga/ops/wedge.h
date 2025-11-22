#pragma once

#include "ga/ops/geometric.h"
#include "ga/multivector.h"

namespace ga::ops {

using ga::Multivector;

inline bool keepWedgeGrade(const int gradeA, const int gradeB, const int gradeR) {
    return gradeR == gradeA + gradeB;
}

// Outer product of two multivectors
inline Multivector wedge(const Multivector& A, const Multivector& B) {
    const Multivector result{ga::ops::geometricProductFiltered(A, B, &keepWedgeGrade)};
    return result;
}

} // namespace ga::ops
