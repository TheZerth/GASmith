GASmith is a light weight Geometric Algebra library for C++.

Intended for use by beginners to GA for tinkeringand learning.

I noticed while trying to teach myself GA as a non-mathematician that most libraries were too heavyweight.
Full of complex math jargon and too much code that made it hard to get started.

So I made this.

Easily add GA to your project. Create and manipulate multivectors with ease.
It's that simple.

Requirements:
- C++17
- CMake
- Python (If you want automated unit tests and benchmarks)

Usage:
1. Clone the repo
2. Run CMake to build the library and examples
3. Include GASmith in your project and start using Geometric Algebra!

Usage Syntax:
```
#include <ga/e3.h>
#include <ga/ops/operators.h>

using namespace ga;
using namespace ga::e3;

int main() {
    auto a = scalar(1.0f) + 2.0f*e1 + 3.0f*e23;
    auto b =          -e2 + 0.5f*e123;

    auto gp = a * b;
    auto out = a ^ b;
    auto inn = a & b;

    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "a*b = " << gp << "\n";
}
```

Core Syntax:
```
#include <GASmith.h>

auto sig = Signature(3,0,0,true);
Algebra alg(sig);

Multivector e1 = ...; // basis e1
Multivector e2 = ...; // basis e2
Rotor R = Rotor::fromPlaneAngle(alg, e1, e2, angle);

Multivector v = ...; // some vector
Multivector v_rot = R.apply(v);
```

Contributing:
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a pull request

Notes:
- the run_benchmarks build target runs benchmarks for all major functions and uploads to influxDb. This allows for automated tracking of performance changes over time.
- the run_tests build target runs unit tests. This will ensure none of your changes break anything.
- benchmarks can also be run through GitHub Actions.
