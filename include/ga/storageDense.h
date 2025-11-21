// C++
#include <vector>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <cstdlib>

namespace ga {

    struct DenseStorage {
        std::vector<double> coefficients;
        uint32_t dimensions;

        explicit DenseStorage(uint32_t dims) : dimensions(dims) {
            if (dims >= (sizeof(size_t) * 8)) {
                std::cerr << "DenseStorage: dims too large: " << dims << "\n";
                std::abort();
            }
            size_t required_size = 1ULL << dims;
            coefficients.resize(required_size, 0.0);
        }

        double& operator[](uint32_t mask) {
            // Diagnostic: print offending values before aborting
            if (static_cast<size_t>(mask) >= coefficients.size()) {
                std::cerr << "DenseStorage::operator[] - mask out of range\n"
                          << "  mask = " << mask << "\n"
                          << "  coefficients.size() = " << coefficients.size() << "\n"
                          << "  dimensions = " << dimensions << "\n";
                std::abort();
            }
            return coefficients[mask];
        }

        const double& operator[](uint32_t mask) const {
            if (static_cast<size_t>(mask) >= coefficients.size()) {
                std::cerr << "DenseStorage::operator[] (const) - mask out of range\n"
                          << "  mask = " << mask << "\n"
                          << "  coefficients.size() = " << coefficients.size() << "\n"
                          << "  dimensions = " << dimensions << "\n";
                std::abort();
            }
            return coefficients[mask];
        }

        size_t size() const {
            return coefficients.size();
        }
    };

}
