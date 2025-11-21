#pragma once
#include "algebra.h"
#include "storageDense.h"

using ga::Algebra;
using ga::DenseStorage;

namespace ga {

struct Multivector {
    const Algebra* alg;     // pointer to algebra descriptor
    DenseStorage storage;   // coefficients indexed by mask

    Multivector(const Algebra& a)
        : alg(&a), storage(a.dimensions) {
    }

    double component(BladeMask m) const { return storage[m]; }
    void setComponent(BladeMask m, double value) { storage[m] = value; }
};

}

