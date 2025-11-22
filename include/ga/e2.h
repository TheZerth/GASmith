#pragma once
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/basis.h"

namespace ga::e2 {

    // Euclidean 2D signature (+,+,+)
    inline const Signature signature{2, 0, 0, true};
    inline const Algebra   algebra{signature};

    // Helpers to construct basis blades in this algebra
    inline Multivector scalar(float s) {
        Multivector mv(algebra);
        mv.setComponent(static_cast<BladeMask>(0), s);
        return mv;
    }

    inline Multivector basis(int axisIndex) {
        Multivector mv(algebra);
        mv.setComponent(Blade::getBasis(axisIndex), 1.0f);
        return mv;
    }

    // Named basis vectors
    inline const Multivector e1 = basis(0);
    inline const Multivector e2 = basis(1);

    // Named bivectors
    inline Multivector bivector(int i, int j) {
        Multivector mv(algebra);
        BladeMask m = static_cast<BladeMask>(Blade::getBasis(i) | Blade::getBasis(j));
        mv.setComponent(m, 1.0f);
        return mv;
    }

    inline const Multivector e12 = bivector(0, 1);

} // namespace ga::e3
