#pragma once

namespace ga::ops {

    using ga::Algebra;
    using ga::Multivector;
    using ga::BladeMask;
    using ga::Blade;

// -------------- Reverse (~A) --------------------------------------------
//
// For each basis blade of grade r: sign = (-1)^{r(r-1)/2}
inline Multivector reverse(const Multivector& A) {
    const Algebra* alg = A.alg;
    if (!alg) {
        // You can throw or assert here if you want
        return A;
    }

    const int dims = alg->dimensions;
    const int bladeCount = 1 << dims;

    Multivector result(*alg);

    for (int i = 0; i < bladeCount; ++i) {
        const auto m = static_cast<BladeMask>(i);
        const double c = A.component(m);
        if (c == 0.0)
            continue;

        const int r = Blade::getGrade(m);
        // r(r-1)/2 mod 2:
        const int exponent = (r * (r - 1) / 2) & 1;
        const double sign = (exponent == 0) ? 1.0 : -1.0;

        result.setComponent(m, static_cast<float>(c * sign));
    }

    return result;
}

// -------------- Grade involution (A^ = sum (-1)^r A_r) ------------------
//
// For each grade-r blade: sign = (-1)^r
inline Multivector gradeInvolution(const Multivector& A) {
    const Algebra* alg = A.alg;
    if (!alg) {
        return A;
    }

    const int dims = alg->dimensions;
    const int bladeCount = 1 << dims;

    Multivector result(*alg);

    for (int i = 0; i < bladeCount; ++i) {
        const auto m = static_cast<BladeMask>(i);
        const double c = A.component(m);
        if (c == 0.0)
            continue;

        const int r = Blade::getGrade(m);
        const double sign = (r & 1) ? -1.0 : 1.0;  // (-1)^r

        result.setComponent(m, static_cast<float>(c * sign));
    }

    return result;
}

// -------------- Clifford conjugation (combination of both) --------------
//
// For each grade-r blade: sign = (-1)^{r(r+1)/2}
inline Multivector cliffordConjugate(const Multivector& A) {
    const Algebra* alg = A.alg;
    if (!alg) {
        return A;
    }

    const int dims = alg->dimensions;
    const int bladeCount = 1 << dims;

    Multivector result(*alg);

    for (int i = 0; i < bladeCount; ++i) {
        const auto m = static_cast<BladeMask>(i);
        const double c = A.component(m);
        if (c == 0.0)
            continue;

        const int r = ga::Blade::getGrade(m);
        const int exponent = (r * (r + 1) / 2) & 1;
        const double sign = (exponent == 0) ? 1.0 : -1.0;

        result.setComponent(m, static_cast<float>(c * sign));
    }

    return result;
}

} // namespace ga::ops
