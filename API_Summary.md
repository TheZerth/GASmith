
This document describes the current architecture and public API of your library as implemented in:

- `signature.h`
- `algebra.h`
- `basis.h`
- `storageDense.h`
- `multivector.h`
- `linearMap.h`
- `versor.h`
- `rotor.h`
- `policies.h`
- `geometric.h`
- `ops/blade.h`
- `ops/wedge.h`
- `ops/inner.h`
- `ops/dual.h`
- `ops/involutions.h`

## 1. Global Config & Policies

### 1.1 `ga::Policies` (`policies.h`)

```cpp
namespace ga {

struct DefaultPolicies {
    using Scalar = float;

    static constexpr Scalar epsilon() noexcept {
        return static_cast<Scalar>(1e-6);
    }
};

using Policies = DefaultPolicies;

} // namespace ga
```

* Defines the **scalar type** used internally (`float`).
* Provides a global **tolerance**, `Policies::epsilon()`, used for:

    * Avoiding exact `== 0` comparisons in normalization (`Versor::inverse`, `Rotor::normalize`).
    * Guarding against inversion/normalization of near-zero norms.

You can later generalize by swapping `Policies` (e.g. `double` or custom eps).

---

## 2. Signature & Algebra Configuration

### 2.1 `ga::Signature` (`signature.h`)

Represents the **metric structure** of a real Clifford algebra (Cl(p, q, r)) with up to 8 axes:

```cpp
namespace ga {

static constexpr int MAX_DIMENSIONS = 8;

using Metric = std::array<int, MAX_DIMENSIONS>;
using Mask   = std::array<bool, MAX_DIMENSIONS>;

struct Signature {
    // counts
    int p_, q_, r_;

    // diagonal metric g_ii ∈ {+1, -1, 0}
    Metric metric_;

    // orientation
    bool isRightHanded_;

    // number of used axes
    int dimensionsUsed_;
    ...
};

} // namespace ga
```

**Semantics:**

* `p()` – number of **positive** axes ((g_{ii} = +1)).
* `q()` – number of **negative** axes ((g_{ii} = -1)).
* `r()` – number of **null** axes ((g_{ii} = 0)), i.e. degeneracy.
* `dimensionsUsed()` – total dimension (n = p+q+r \le MAX_DIMENSIONS).
* `metric()` – the diagonal metric array.
* `isRightHanded()` / `isLeftHanded()` – orientation of the pseudoscalar.

**Metric access:**

```cpp
int  metricLookup(int i, int j) const; // diag metric, off-diagonal = 0, out-of-range = -2
int  getSign(int i) const;             // g_ii
bool isPos(int i) const;
bool isNeg(int i) const;
bool isZero(int i) const;
bool isDegenerate() const;             // (r_ > 0)
```

**Construction modes:**

1. From `(p, q, r, handedness)`:

    * Builds diagonal metric: `+1` for first `p`, `-1` for next `q`, `0` for next `r`.
2. From masks `Mask pMask, qMask, rMask`:

    * Assigns specific indices as +/−/0.
    * Ensures masks don’t overlap.
3. From a raw `Metric` and `axisCount`:

    * Counts `+1` → `p_`, `-1` → `q_`, “anything else” → `r_`.

This gives a flexible but explicit declaration of the Clifford algebra’s metric, including degeneracy.

---

### 2.2 `ga::Algebra` (`algebra.h`)

A lightweight descriptor of a specific algebra instance:

```cpp
namespace ga {

struct Algebra {
    Signature signature;
    int       dimensions;

    Algebra();                         // default Signature, uses signature.dimensionsUsed()
    explicit Algebra(const Signature& sig);

    void setSignature(const Signature& sig);
};

} // namespace ga
```

**Invariant:**

```cpp
dimensions == signature.dimensionsUsed();
```

This type serves as a **shared context** for:

* Multivectors
* Linear maps
* Versors / rotors

and ensures all GA objects agree on the underlying metric and dimension.

---

## 3. Basis & Blade Representation

### 3.1 `ga::BladeMask` & `ga::Blade` (`basis.h`)

```cpp
namespace ga {

static constexpr int MAX_DIMENSIONS = 8;
using BladeMask = std::uint8_t;

struct Blade {
    BladeMask mask{}; // which axes are present
    int       sign{}; // 0 = zero blade, ±1 = oriented blade

    constexpr Blade() = default;
    constexpr Blade(BladeMask m, int s) : mask(m), sign(s) {}
    ...
};

} // namespace ga
```

**Concepts:**

* `mask` – bit `i` set ⇔ basis vector `e_i` present.
* Grade = `popcount(mask)`.
* `sign` – orientation of the blade; `sign == 0` means **zero blade**.
* Scalar basis: `mask == 0 && sign != 0` (unit scalar basis `1`).
* Zero blade: `sign == 0`.

**Key helpers:**

```cpp
int  getGrade(BladeMask mask);
bool hasAxis(BladeMask mask, int i);
BladeMask getBasis(int axisIndex); // mask for e_i

int  highestAxis(BladeMask mask);
bool doesOverlap(BladeMask a, BladeMask b);
bool isZero(Blade blade);
bool isScalarBasis(Blade blade);

BladeMask addAxis(BladeMask mask, int axisIndex);
BladeMask removeAxis(BladeMask mask, int axisIndex);
BladeMask toggleAxis(BladeMask mask, int axisIndex);
```

### 3.2 Canonical blade construction: `Blade::makeBlade`

```cpp
static constexpr Blade makeBlade(const int* basis, int numBasis);
```

* Takes a list of axis indices, sorts them, tracks permutation parity.
* If any index repeats, returns zero blade (`e_i ∧ e_i = 0`).
* Computes:

    * `mask` = OR of `1 << idx` over all axes.
    * `sign` = +1 or −1 based on permutation parity.

This encodes the **exterior basis** of the GA, independent of metric.

### 3.3 Exterior product of basis blades: `Blade::combineBlade`

```cpp
static constexpr Blade combineBlade(Blade a, Blade b);
```

Implements the **wedge product** of basis blades:

* Zero annihilates.
* Scalar basis is identity.
* If `a` and `b` share any axis → wedge is zero.
* Otherwise:

    * Result mask = `a.mask ^ b.mask` (symmetric difference).
    * Sign = parity of the swaps needed to bring the combined basis into ascending index order.

This is the standard exterior algebra rule.

---

## 4. Dense Multivector Storage

### 4.1 `ga::DenseStorage` (`storageDense.h`)

```cpp
namespace ga {

struct DenseStorage {
    static constexpr std::size_t MAX_ELEMENTS = 256; // 2^8

    float   coefficients[MAX_ELEMENTS]{};
    uint8_t dimensions; // number of axes (n ≤ 8)

    explicit DenseStorage(uint8_t dims);

    float&       operator[](std::size_t mask);
    float const& operator[](std::size_t mask) const;

    std::size_t size() const;     // = 1 << dimensions
    static constexpr std::size_t capacity(); // = MAX_ELEMENTS
};

} // namespace ga
```

* Indexing is by `BladeMask` (0..(1<<n)-1).
* `size()` = total number of basis blades for dimension `n`.
* Only indices `< size()` are valid; asserts in debug builds enforce this.

This is a canonical dense representation of a multivector’s coordinates.

---

## 5. Multivector Type

### 5.1 `ga::Multivector` (`multivector.h`)

```cpp
namespace ga {

struct Multivector {
    const Algebra* alg;     // algebra descriptor
    DenseStorage   storage; // coefficients indexed by mask

    explicit Multivector(const Algebra& a)
        : alg(&a), storage(a.dimensions) {}

    double component(BladeMask m) const { return storage[m]; }
    void   setComponent(BladeMask m, double value) { storage[m] = static_cast<float>(value); }
};

} // namespace ga
```

**Semantics:**

* Represents a general GA element:
  [
  A = \sum_{m} a_m E_m,
  ]
  where `m` is a bitmask and `E_m` is the basis blade with that mask.
* `alg` must remain valid for the lifetime of the `Multivector`.
* The algebra pointer ensures you never multiply elements of different algebras (all ops check `A.alg == B.alg`).

---

## 6. Basis-Level Clifford Product

### 6.1 `ga::ops::geometricProductBlade` (`blade.h`)

```cpp
namespace ga::ops {

Blade geometricProductBlade(const Blade& a,
                            const Blade& b,
                            const Signature& sig);

} // namespace ga::ops
```

Implements the **Clifford product** on basis blades (e_A e_B):

1. **Special cases:**

    * Zero annihilates: `0 * anything = 0`.
    * Scalars as identity:

        * `1 * B = B`, `B * 1 = B`, sign handled by `a.sign * b.sign`.

2. **Swap parity (exterior sign):**

    * Counts how many “inversions” appear when concatenating axes from `a` and `b` and sorting into canonical order.
    * Contributes factor ((-1)^{\text{swaps}}).

3. **Metric contraction:**

    * Finds overlapping axes in `am & bm`:

        * For each overlap index `i`, multiplies sign by `sig.getSign(i)`:
          [
          e_i e_i = g_{ii} \in {+1, -1, 0}.
          ]
        * If any `g_{ii} == 0` (null axis) → result = zero blade.

4. **Result mask & sign:**

    * Result mask = `am ^ bm` (symmetric difference).
    * If result mask is zero → scalar blade `{0, sign}`.
    * Otherwise → blade `{resultMask, sign}`.

This is the standard algebra rule:

[
e_A e_B = \left(a.\text{sign} \cdot b.\text{sign} \cdot (-1)^{\text{swaps}}
\cdot \prod_{i \in A\cap B} g_{ii}\right)
, e_{A\Delta B}.
]

---

## 7. Multivector-Level Clifford Product

### 7.1 Filtered product: `ga::ops::geometricProductFiltered` (`geometric.h`)

```cpp
namespace ga::ops {

using GradeFilterFn = bool (*)(int gradeA, int gradeB, int gradeR);

Multivector geometricProductFiltered(const Multivector& A,
                                     const Multivector& B,
                                     GradeFilterFn keep);

Multivector geometricProduct(const Multivector& A,
                             const Multivector& B);

} // namespace ga::ops
```

**Process:**

* Validates that `A.alg` and `B.alg` are non-null and equal.

* For each basis mask `i` and `j`:

    * `coeffA = A.component(i)`, `coeffB = B.component(j)`; skip if zero.
    * Compute `gp = geometricProductBlade(Blade{i, +1}, Blade{j, +1}, A.alg->signature)`.
    * Optionally compute grades `gradeA`, `gradeB`, `gradeR`.
    * If `keep` is provided and returns `false` → skip term.
    * Otherwise, accumulate:
      [
      \text{result}_{gp.mask} \mathrel{+}= coeffA \cdot coeffB \cdot gp.sign.
      ]

* `geometricProduct(A,B)` calls `geometricProductFiltered(A,B,nullptr)` (no grade filter).

This is the pure linear extension of the basis-level Clifford product.

---

## 8. Derived Products

### 8.1 Outer (wedge) product (`wedge.h`)

```cpp
namespace ga::ops {

bool keepWedgeGrade(int gradeA, int gradeB, int gradeR) {
    return gradeR == gradeA + gradeB;
}

inline Multivector wedge(const Multivector& A, const Multivector& B) {
    return geometricProductFiltered(A, B, &keepWedgeGrade);
}

} // namespace ga::ops
```

Implements
[
A \wedge B = \sum_{r,s} \langle A_r B_s \rangle_{r+s},
]
in line with modern GA.

### 8.2 Inner product & contractions (`inner.h`)

```cpp
namespace ga::ops {

// Hestenes inner: gradeR == |gradeA - gradeB|
bool keepInnerGrade(int gradeA, int gradeB, int gradeR);

Multivector inner(const Multivector& A, const Multivector& B);

// Left contraction: gradeR == (gradeB - gradeA) if gradeA <= gradeB
bool keepLeftContractionGrade(int gradeA, int gradeB, int gradeR);
Multivector leftContraction(const Multivector& A, const Multivector& B);

// Right contraction: gradeR == (gradeA - gradeB) if gradeA >= gradeB
bool keepRightContractionGrade(int gradeA, int gradeB, int gradeR);
Multivector rightContraction(const Multivector& A, const Multivector& B);

} // namespace ga::ops
```

These match standard GA definitions based on grade projection of the Clifford product.

---

## 9. Involutions

### 9.1 Reversion (`reverse`) (`involutions.h`)

```cpp
namespace ga::ops {

Multivector reverse(const Multivector& A);

} // namespace ga::ops
```

* For each blade of grade `r`, applies sign:
  [
  (-1)^{r(r-1)/2}.
  ]

### 9.2 Grade involution (`gradeInvolution`)

```cpp
Multivector gradeInvolution(const Multivector& A);
```

* For grade `r`, sign:
  [
  (-1)^r.
  ]

### 9.3 Clifford conjugation (`cliffordConjugate`)

```cpp
Multivector cliffordConjugate(const Multivector& A);
```

* Combined effect:
  [
  (-1)^{r(r+1)/2},
  ]
  i.e. reverse followed by grade involution (or vice versa).

All three are canonical in the GA literature.

---

## 10. Dual / Hodge Dual

### 10.1 `ga::ops::dual` (`dual.h`)

```cpp
namespace ga::ops {

Multivector dual(const Multivector& A);

} // namespace ga::ops
```

Algorithm:

1. Let `dims = alg->dimensions`. Define pseudoscalar mask:

   ```cpp
   BladeMask I_mask = static_cast<BladeMask>((1u << dims) - 1u);
   ```
2. For each non-zero component of `A` with mask `m`:

    * Complement mask:

      ```cpp
      BladeMask comp = I_mask ^ m;
      ```
    * Compute `gp = geometricProductBlade(Blade{m, +1}, Blade{comp, +1}, signature)`.
    * **Updated behavior:**

        * If `gp.sign == 0` **or** `gp.mask != I_mask`, treat the dual as undefined for this blade and **skip** it (degenerate or ill-defined in this metric).
        * Otherwise, contribute:
          [
          \text{dual}(A)_{\text{comp}} \mathrel{+}= a_m \cdot gp.sign.
          ]

This yields the standard dual in non-degenerate algebras and clearly surfaces degeneracy instead of silently misbehaving.

---

## 11. Linear Maps & Outermorphisms

### 11.1 `ga::LinearMap` (`linearMap.h`)

```cpp
namespace ga {

struct LinearMap {
    const Algebra* alg = nullptr;
    float          m[8][8]; // L(e_col) = sum_row m[row][col] e_row

    LinearMap();                          // zero map, no algebra
    explicit LinearMap(const Algebra& a); // identity on a

    static LinearMap identity(const Algebra& a);
    static LinearMap zero(const Algebra& a);

    void  set(int row, int col, float value);
    float get(int row, int col) const;

    Multivector applyToVector(const Multivector& v) const;
    Multivector apply(const Multivector& A) const;
};

} // namespace ga
```

**`applyToVector(v)`**:

* Treats `v` as a pure **grade-1** multivector.
* Extracts `v_j` from basis vectors `e_j`.
* Computes `w_i = ∑_j m[i][j] v_j`.
* Returns a new vector multivector.

**`apply(A)` (outermorphism):**

1. Precompute `vecImages[j] = L(e_j)` via `applyToVector`.
2. Precompute `bladeImages[mask]` for all basis blades:

    * `mask == 0` → scalar 1: `L(1) = 1`.
    * grade-1: `L(e_j) = vecImages[j]`.
    * higher grade:

        * Decompose mask as `firstAxis` + `remaining`.
        * Recursively: `L(mask) = L(e_firstAxis) ∧ L(remaining)`.
3. Extend linearly:
   [
   L(A) = \sum_m a_m L(E_m).
   ]

This is the standard **outermorphism** extension of a linear map (L: V \to V) to the full exterior/Clifford algebra.

---

## 12. Versors

### 12.1 `ga::Versor` (`versor.h`)

```cpp
namespace ga {

class Versor {
public:
    Multivector mv;

    Versor() = default;
    explicit Versor(const Multivector& v);

    Versor(const Algebra& algebra, const Multivector& v);

    const Algebra* algebra() const noexcept;
    bool           isValid() const noexcept;

    Multivector inverse() const;
    Multivector apply(const Multivector& X) const;
};

} // namespace ga
```

**Semantics:**

* A `Versor` is an **invertible multivector** that acts on other multivectors by:

  [
  X' = V X V^{-1}.
  ]

* Inverse is defined via reversion:

  ```cpp
  Multivector vrev = reverse(mv);
  Multivector norm2_mv = geometricProduct(mv, vrev);
  float s = (scalar part of norm2_mv);
  ```

  **Updated behavior:**

    * Use `Policies::epsilon()`:

      ```cpp
      if (std::fabs(s) <= Policies::epsilon()) {
          throw std::runtime_error("scalar norm (V ~V) is too close to zero");
      }
      ```
    * Then:
      [
      V^{-1} = \frac{\tilde{V}}{s}.
      ]

* `apply(X)`:

    * Ensures `mv.alg` and `X.alg` are non-null and equal.
    * Computes:

      ```cpp
      Multivector invV = inverse();
      return geometricProduct( geometricProduct(mv, X), invV );
      ```

So versors are built using the canonical GA formula for invertible multivectors.

---

## 13. Rotors (Even Versors for Rotations/Boosts)

### 13.1 `ga::Rotor` (`rotor.h`)

```cpp
namespace ga {

class Rotor {
public:
    Multivector mv;

    Rotor() = default;
    explicit Rotor(const Multivector& R);

    const Algebra*   algebra() const noexcept;
    bool             isValid() const noexcept;
    Multivector&     value() noexcept;
    const Multivector& value() const noexcept;

    void        normalize();
    Multivector apply(const Multivector& X) const;

    static Rotor fromBivectorAngle(const Multivector& B, float theta);
    static Rotor fromPlaneAngle(const Multivector& a,
                                const Multivector& b,
                                float theta);
};

} // namespace ga
```

**`normalize()`**:

* Computes `R ~R` and extracts scalar part `s`.
* **Updated behavior:**

    * Uses `Policies::epsilon()`:

      ```cpp
      if (std::fabs(s) <= Policies::epsilon()) {
          throw std::runtime_error("rotor norm^2 is too close to zero");
      }
      ```

    * Scales all coefficients by (1 / \sqrt{|s|}) so that (R\tilde{R} ≈ 1).

**`apply(X)`**:

* Checks algebra compatibility.
* Uses the **unit rotor sandwich**:

  [
  X' = R X \tilde{R}.
  ]

**`fromBivectorAngle(B, θ)`**:

* Assumes `B` is (approximately) a unit bivector representing a plane.

* Builds:

  [
  R = \cos(\theta/2) - B \sin(\theta/2),
  ]

* Then calls `normalize()` for safety.

**`fromPlaneAngle(a, b, θ)`**:

* Constructs plane bivector:

  ```cpp
  Multivector B = wedge(a, b);
  ```

* **Updated metric-aware normalization:**

  ```cpp
  Multivector BB   = inner(B, B); // B ⋅ B
  float norm2      = scalar_part(BB);
  auto  eps        = Policies::epsilon();

  if (std::fabs(norm2) <= eps) {
      throw std::runtime_error("a ∧ b has zero (or near-zero) norm");
  }

  float inv_mag = 1.0f / std::sqrt(std::fabs(norm2));

  for (i: 0..N-1)
      B.storage[i] *= inv_mag;
  ```

* Then uses `fromBivectorAngle(B, θ)`.

This ensures rotor construction respects the **metric structure** of the algebra instead of relying on a naive Euclidean coefficient norm.

---

## 14. Axioms & Design Guarantees

1. **Clifford product is explicit and standard:**

    * Basis-level rule:
      [
      e_A e_B = (-1)^{\text{swaps}} \prod_{i\in A\cap B} g_{ii} , e_{A\Delta B}
      ]
      with null directions correctly annihilating products.
    * Multivector-level product is the linear extension of this rule.

2. **Exterior, inner, and contraction products are grade projections:**

    * No bespoke or ad-hoc formulas — they’re all defined in terms of grade-filtered `geometricProduct`.

3. **Involutions are textbook:**

    * Reverse: ((-1)^{r(r-1)/2}).
    * Grade involution: ((-1)^r).
    * Clifford conjugation: ((-1)^{r(r+1)/2}).

4. **Dual is defined via complement + Clifford sign:**

    * Correct for non-degenerate metrics.
    * Degenerate signatures are handled explicitly by skipping ill-defined cases.

5. **Linear maps extend via outermorphism:**

    * (L(a \wedge b) = L(a) \wedge L(b)).
    * Scalars fixed: (L(1) = 1).

6. **Versors and rotors use canonical GA formulas:**

    * Versor inverse: (\tilde{V} / (V \tilde{V})_{\text{scalar}}).
    * Rotor action: (X' = R X \tilde{R}).
    * Rotor construction from bivectors/planes follows standard (R = \cos(\theta/2) - B \sin(\theta/2)), with proper metric-aware normalization.

7. **Numeric robustness is centralized via `Policies::epsilon()`:**

    * Inversion and normalization operations throw when norms are **too close to zero**, not only exactly zero.

---