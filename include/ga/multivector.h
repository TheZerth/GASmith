#pragma once
#include "storageDense.h"
#include "algebra.h"
#include "basis.h"
// A mutlivector is our generalized maths object. It essentially is a combination of a BladeMask and an Algebra metric.
// The BladeMask defines what axis the multivector contains.
// The coeffiecients for this mask are coded in the denseStorage which is an array of 256 ints.
// This supports a maximum of 8D maths for any clifford metric.
// The algebra signature is needed to ensure that Multivectors from different metrics aren't computed together.
// Aditionally the algebra metric helps define how axis overlaps behave.
// An algebra pointer is used as the intention is to pre-compute algebra lookup tables and all multivectors of that metric will have access to it.

namespace ga {

    using ga::Algebra;
    using ga::DenseStorage;
    using ga::BladeMask;

    struct Multivector {
        const Algebra* alg;     // pointer to algebra descriptor
        DenseStorage storage;   // coefficients indexed by mask

        explicit Multivector(const Algebra& a)
            : alg(&a), storage(a.dimensions) {}

        [[nodiscard]] double component(const BladeMask m) const { return storage[m]; }
        void setComponent(const BladeMask m, const double value) { storage[m] = value; }

    };

    inline Multivector operator+(const Multivector& A, const Multivector& B) {
        if (!A.alg || !B.alg || A.alg != B.alg) {
            throw std::invalid_argument("ga::operator+: operands must share the same Algebra");
        }
        //Grab current algebra
        Multivector result(*A.alg);
        //Grab dimensions
        const int dims = A.alg->dimensions;
        //Determine number of coefficients
        const std::size_t n = (1u << dims);
        //Iterate over all blades and add them
        for (int i = 0; i < n; ++i) {
            result.storage[i] = A.storage[i] + B.storage[i];
        }
        return result;
    }

    inline Multivector operator-(const Multivector& A, const Multivector& B) {
        if (!A.alg || !B.alg || A.alg != B.alg) {
            throw std::invalid_argument("ga::operator+: operands must share the same Algebra");
        }
        //Grab current algebra
        Multivector result(*A.alg);
        //Grab dimensions
        const int dims = A.alg->dimensions;
        //Determine number of coefficients
        const std::size_t n = (1u << dims);
        //Iterate over all blades and add them
        for (int i = 0; i < n; ++i) {
            result.storage[i] = A.storage[i] - B.storage[i];
        }
        return result;
    }

    inline Multivector operator*(const float f, const Multivector& A) {
        if (!A.alg) {
            throw std::invalid_argument("ga::operator+: multivector has no Algebra");
        }
        //Grab current algebra
        Multivector result(*A.alg);
        //Grab dimensions
        const int dims = A.alg->dimensions;
        //Determine number of coefficients
        const std::size_t n = (1u << dims);
        //Iterate over all blades and multiply them
        for (int i = 0; i < n; ++i) {
            result.storage[i] = A.storage[i] * f;
        }
        return result;
    }

    inline Multivector operator*(const Multivector& A, const float f) {
        return f * A;
    }

}

