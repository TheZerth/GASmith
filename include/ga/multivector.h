#pragma once
#include "algebra.h"
#include "storageDense.h"

// A mutlivector is our generalized maths object. It essentially is a combination of a BladeMask and an Algebra metric.
// The BladeMask defines what axis the multivector contains.
// The coeffiecients for this mask are coded in the denseStorage which is an array of 256 ints.
// This supports a maximum of 8D maths for any clifford metric.
// The algebra signature is needed to ensure that Multivectors from different metrics aren't computed together.
// Aditionally the algebra metric helps define how axis overlaps behave.
// An algebra pointer is used as the intention is to pre-compute algebra lookup tables and all multivectors of that metric will have access to it.

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

