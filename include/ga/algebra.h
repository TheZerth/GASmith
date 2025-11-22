// C++
#pragma once
#include "ga/signature.h"

namespace ga {

    using ga::Signature;

    struct Algebra {
        Signature signature;
        int dimensions;

        // default ctor initializes dimensions from signature
        Algebra() : signature(), dimensions(signature.dimensionsUsed()) {}

        // construct with a Signature and sync dimensions
        explicit Algebra(const Signature& sig) : signature(sig), dimensions(sig.dimensionsUsed()) {}

        // helper to update signature and keep dimensions in sync
        void setSignature(const Signature& sig) {
            signature = sig;
            dimensions = sig.dimensionsUsed();
        }
    };

} // namespace ga
