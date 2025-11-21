#pragma once

Multivector geometricProduct(const Multivector& A, const Multivector& B) {
    Multivector out(*A.alg);
    for (BladeMask i = 0; i < (1u << A.alg->dimensions); ++i) {
        double a = A.storage[i];
        if (a == 0) continue;

        Blade ba{i, +1};

        for (BladeMask j = 0; j < (1u << A.alg->dimensions); ++j) {
            double b = B.storage[j];
            if (b == 0) continue;

            Blade bb{j, +1};
            Blade bc = geometricProductBlade(ba, bb, A.alg->signature);

            if (!Blade::isZero(bc)) {
                out.storage[bc.mask] += a * b * bc.sign;
            }
        }
    }
    return out;
}