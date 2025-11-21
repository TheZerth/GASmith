#pragma once

#include <bit>        // std::popcount
#include <stdexcept>  // std::invalid_argument
#include <cstdlib>    // std::abs

#include "ga/multivector.h"
#include "ga/algebra.h"
#include "ga/signature.h"
#include "ga/basis.h"
#include "ga/ops/blade.h"
#include "ga/ops/geometric.h"

namespace ga::ops {

    // ---------------- Hestenes inner product filter ---------------------------
    // keep gradeR == |gradeA - gradeB|
    inline bool keepInnerGrade(int gradeA, int gradeB, int gradeR) {
        int diff = std::abs(gradeA - gradeB);
        return gradeR == diff;
    }

    // ---------------- Left contraction filter ---------------------------------
    // keep gradeR == (gradeB - gradeA) if gradeA <= gradeB
    inline bool keepLeftContractionGrade(int gradeA, int gradeB, int gradeR) {
        if (gradeA > gradeB)
            return false;
        int expected = gradeB - gradeA;
        return gradeR == expected;
    }

    // ---------------- Right contraction filter --------------------------------
    // keep gradeR == (gradeA - gradeB) if gradeA >= gradeB
    inline bool keepRightContractionGrade(int gradeA, int gradeB, int gradeR) {
        if (gradeA < gradeB)
            return false;
        int expected = gradeA - gradeB;
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