#pragma once
#include "ga/signature.h"
#include "ga/algebra.h"
#include "ga/multivector.h"
#include "ga/basis.h"

namespace ga::e3 {

    // Euclidean 3D signature (+,+,+)
    inline const Signature signature{3, 0, 0, true};
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
    inline const Multivector e3 = basis(2);

    // Named bivectors
    inline Multivector bivector(int i, int j) {
        Multivector mv(algebra);
        BladeMask m = static_cast<BladeMask>(Blade::getBasis(i) | Blade::getBasis(j));
        mv.setComponent(m, 1.0f);
        return mv;
    }

    inline const Multivector e12 = bivector(0, 1);
    inline const Multivector e13 = bivector(0, 2);
    inline const Multivector e23 = bivector(1, 2);

    // Trivector
    inline const Multivector e123 = []{
        Multivector mv(algebra);
        BladeMask m = static_cast<BladeMask>(
            Blade::getBasis(0) | Blade::getBasis(1) | Blade::getBasis(2));
        mv.setComponent(m, 1.0f);
        return mv;
    }();

} // namespace ga::e3
