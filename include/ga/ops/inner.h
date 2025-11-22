#pragma once

#include <cstdlib>    // std::abs
#include "ga/multivector.h"

namespace ga::ops {

    using ga::Multivector;

    // ---------------- Hestenes inner product filter ---------------------------
    // keep gradeR == |gradeA - gradeB|
    inline bool keepInnerGrade(const int gradeA, const int gradeB, const int gradeR) {
        const int diff = std::abs(gradeA - gradeB);
        return gradeR == diff;
    }

    // ---------------- Left contraction filter ---------------------------------
    // keep gradeR == (gradeB - gradeA) if gradeA <= gradeB
    inline bool keepLeftContractionGrade(const int gradeA, const int gradeB, const int gradeR) {
        if (gradeA > gradeB)
            return false;
        const int expected = gradeB - gradeA;
        return gradeR == expected;
    }

    // ---------------- Right contraction filter --------------------------------
    // keep gradeR == (gradeA - gradeB) if gradeA >= gradeB
    inline bool keepRightContractionGrade(const int gradeA, const int gradeB, const int gradeR) {
        if (gradeA < gradeB)
            return false;
        const int expected = gradeA - gradeB;
        return gradeR == expected;
    }

    // --------------------------------------------------------------------------
    //  Hestenes inner product: A ⋅ B
    // --------------------------------------------------------------------------
    inline Multivector inner(const Multivector& A, const Multivector& B) {
        return ga::ops::geometricProductFiltered(A, B, &keepInnerGrade);
    }

    // --------------------------------------------------------------------------
    //  Left contraction: A ⟍ B
    // --------------------------------------------------------------------------
    inline Multivector leftContraction(const Multivector& A, const Multivector& B) {
        return ga::ops::geometricProductFiltered(A, B, &keepLeftContractionGrade);
    }

    // --------------------------------------------------------------------------
    //  Right contraction: A ⟎ B
    // --------------------------------------------------------------------------
    inline Multivector rightContraction(const Multivector& A, const Multivector& B) {
        return ga::ops::geometricProductFiltered(A, B, &keepRightContractionGrade);
    }

} // namespace ga::ops