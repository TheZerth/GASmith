#pragma once

namespace ga {

struct DenseStorage {
    std::vector<double> data;

    DenseStorage(int dims)
        : data(1 << dims, 0.0) {}

    double& operator[](size_t i) { return data[i]; }
    const double& operator[](size_t i) const { return data[i]; }

    void clear() {
        std::fill(data.begin(), data.end(), 0.0);
    }
};

}
