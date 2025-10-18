cmake_minimum_required(VERSION 4.0)
project(GASmith)

set(CMAKE_CXX_STANDARD 20)

add_library(GASmith STATIC library.cpp
        include/ga/signature.h
        include/ga/policies.h
        include/ga/algebra.h
        include/ga/basis.h
        include/ga/storageDense.h
        include/ga/multivector.h
        include/ga/ops/geometric.h
        include/ga/ops/wedge.h
        include/ga/ops/inner.h
        include/ga/ops/involutions.h
        include/ga/ops/dual.h
        include/ga/linearMap.h
        include/ga/versor.h
        include/ga/rotor.h)
