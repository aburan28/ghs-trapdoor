/-
  GhsTrapdoor.Basic

  A small, mathlib-free formalization of the algebra underlying the GHS
  trap-door system (Gaudry–Hess–Smart, "Constructive and Destructive Facets
  of Weil Descent on Elliptic Curves", J. Cryptology 15(1), 2002):

  * `asp` – the Artin–Schreier operator ℘(t) = t² + t over a commutative ring
    of characteristic 2.  It is additive (`asp_add`), its kernel on a domain
    is exactly {0, 1} (`asp_eq_zero`), and its fibres are the pairs
    {t, t+1} (`asp_fiber_shift`).  In the GHS construction these facts give
    the twist ε ∈ {0,1} selecting the two irreducible components of the Weil
    restriction, and the hyperelliptic involution ŝ ↦ ŝ + 1 of the
    descended curve.

  * `frobSum`/`trace_asp_zero` – the trace identity
      Σ_{i<n} ℘(t)^{2^i} = t^{2^n} + t = 0
    over a field with 2^n elements: the trace-zero condition that makes the
    Artin–Schreier solver in the Sage implementation well defined.

  * `tau_on` – the characteristic-2 curve involution (x,y) ↦ (x, y+x) of
    E : y² + xy = x³ + β, which generates the two components of the Weil
    restriction and the hyperelliptic structure of the descent curve.
-/

universe u

/-- A commutative ring with one (axiomatised from scratch, no mathlib). -/
class CRing (R : Type u) where
  add : R → R → R
  mul : R → R → R
  zero : R
  one : R
  add_assoc (a b c : R) : add (add a b) c = add a (add b c)
  add_comm (a b : R) : add a b = add b a
  add_zero (a : R) : add a zero = a
  mul_assoc (a b c : R) : mul (mul a b) c = mul a (mul b c)
  mul_comm (a b : R) : mul a b = mul b a
  mul_one (a : R) : mul a one = a
  left_distrib (a b c : R) : mul a (add b c) = add (mul a b) (mul a c)

instance {R : Type u} [CRing R] : Add R := ⟨CRing.add⟩
instance {R : Type u} [CRing R] : Mul R := ⟨CRing.mul⟩
instance {R : Type u} [CRing R] : Zero R := ⟨CRing.zero⟩
instance {R : Type u} [CRing R] : One R := ⟨CRing.one⟩

/-- A commutative ring of characteristic two. -/
class CharTwoRing (R : Type u) extends CRing R where
  add_self (a : R) : a + a = 0

section
variable {R : Type u} [CharTwoRing R]

/-! ### Ring plumbing (notation form) -/

theorem add_assoc (a b c : R) : a + b + c = a + (b + c) := CRing.add_assoc a b c
theorem add_comm (a b : R) : a + b = b + a := CRing.add_comm a b
theorem mul_assoc (a b c : R) : a * b * c = a * (b * c) := CRing.mul_assoc a b c
theorem mul_comm (a b : R) : a * b = b * a := CRing.mul_comm a b
theorem mul_one' (a : R) : a * 1 = a := CRing.mul_one a
theorem add_zero' (a : R) : a + 0 = a := CRing.add_zero a
theorem left_distrib (a b c : R) : a * (b + c) = a * b + a * c :=
  CRing.left_distrib a b c

theorem add_left_comm (a b c : R) : a + (b + c) = b + (a + c) := by
  rw [← add_assoc, add_comm a b, add_assoc]

@[simp] theorem add_self (a : R) : a + a = 0 := CharTwoRing.add_self a

@[simp] theorem zero_add (a : R) : 0 + a = a := by
  rw [add_comm, add_zero']

@[simp] theorem one_mul (a : R) : 1 * a = a := by
  rw [← mul_one' a, mul_comm, mul_one']

theorem right_distrib (a b c : R) : (b + c) * a = b * a + c * a := by
  rw [mul_comm (b + c) a, left_distrib, mul_comm a b, mul_comm a c]

@[simp] theorem zero_mul (a : R) : 0 * a = 0 := by
  have h1 : (0 + 0) * a = 0 * a + 0 * a := right_distrib a 0 0
  have h2 : (0 + 0) * a = 0 * a := by rw [add_self]
  rw [← h2, h1, add_self]

/-! ### The Artin–Schreier operator ℘(t) = t² + t -/

/-- The Artin–Schreier operator. -/
def asp (t : R) : R := t * t + t

/-- Frobenius squares are additive in characteristic two. -/
@[simp] theorem sq_add (a b : R) : (a + b) * (a + b) = a * a + b * b := by
  have h1 : (a + b) * (a + b) = (a + b) * a + (a + b) * b := left_distrib (a + b) a b
  have h2 : (a + b) * a = a * a + b * a := right_distrib a a b
  have h3 : (a + b) * b = a * b + b * b := right_distrib b a b
  rw [h1, h2, h3, add_assoc, mul_comm b a, ← add_assoc (a * b) (a * b) (b * b),
    add_self, zero_add]

/-- ℘ is additive: the Frobenius is a ring homomorphism in characteristic 2. -/
theorem asp_add (a b : R) : asp (a + b) = asp a + asp b := by
  simp only [asp, sq_add]
  simp [add_assoc, add_comm, add_left_comm]

@[simp] theorem asp_zero : asp (0 : R) = 0 := by
  have h : asp 0 = (0 : R) * 0 + 0 := rfl
  rw [h, zero_mul, add_self]

@[simp] theorem asp_one : asp (1 : R) = 0 := by
  have h : asp 1 = (1 : R) * 1 + 1 := rfl
  rw [h, one_mul, add_self]

/-- Fibres of ℘ are the pairs {t, t+1}: the twist ε of the GHS components. -/
theorem asp_fiber_shift (c t : R) (h : asp t = c) : asp (t + 1) = c := by
  rw [asp_add, h, asp_one, add_zero']

/-- A characteristic-2 domain. -/
class CharTwoDomain (R : Type u) extends CharTwoRing R where
  mul_eq_zero (a b : R) : a * b = 0 → a = 0 ∨ b = 0

variable {D : Type u} [CharTwoDomain D]

/-- The kernel of ℘ on a characteristic-2 domain is exactly {0, 1}. -/
theorem asp_eq_zero (t : D) (h : asp t = 0) : t = 0 ∨ t = 1 := by
  have hfac : t * (t + 1) = 0 := by
    rw [left_distrib t t 1, mul_one', ← asp, h]
  rcases CharTwoDomain.mul_eq_zero t (t + 1) hfac with h1 | h2
  · exact Or.inl h1
  · refine Or.inr ?_
    have h5 : t + 1 = 1 + 1 := by rw [h2, add_self]
    have h6 : (t + 1) + 1 = (1 + 1) + 1 := by rw [h5]
    rw [add_assoc t 1 1, add_self (1 : D), add_zero', zero_add] at h6
    exact h6

/-! ### Frobenius powers and the trace identity -/

/-- `x ^ n` by iterated multiplication. -/
def npow {R : Type u} [CRing R] (x : R) : Nat → R
  | 0 => 1
  | n + 1 => npow x n * x

instance {R : Type u} [CRing R] : HPow R Nat R := ⟨npow⟩

@[simp] theorem npow_zero (x : R) : x ^ 0 = 1 := rfl

theorem npow_succ (x : R) (n : Nat) : x ^ (n + 1) = x ^ n * x := by
  have h : npow x (n + 1) = npow x n * x := by rw [npow]
  exact h

@[simp] theorem npow_one (x : R) : x ^ 1 = x := by
  rw [npow_succ, npow_zero, one_mul]

theorem pow_add (x : R) : ∀ m n : Nat, x ^ (m + n) = x ^ m * x ^ n := by
  intro m n
  induction n with
  | zero => rw [Nat.add_zero, npow_zero, mul_one']
  | succ n ih =>
      rw [← Nat.add_assoc, Nat.add_comm m n, npow_succ, Nat.add_comm n m, ih,
        npow_succ, mul_assoc]

/-- `(x * y) ^ 2^n = x ^ 2^n * y ^ 2^n`. -/
theorem pow_two_pow_mul (x y : R) : ∀ n : Nat, (x * y) ^ (2 ^ n) = x ^ (2 ^ n) * y ^ (2 ^ n) := by
  intro n
  induction n with
  | zero => simp [Nat.pow_zero, npow_one]
  | succ n ih =>
      have h2 : 2 ^ (n + 1) = 2 ^ n + 2 ^ n := by rw [Nat.pow_succ, Nat.mul_two]
      rw [h2, pow_add, pow_add, pow_add, ih]
      rw [mul_assoc, mul_comm (y ^ 2 ^ n) (x ^ 2 ^ n * y ^ 2 ^ n),
        mul_assoc (x ^ 2 ^ n) (y ^ 2 ^ n) (y ^ 2 ^ n),
        mul_assoc (x ^ 2 ^ n) (x ^ 2 ^ n) (y ^ 2 ^ n * y ^ 2 ^ n)]

/-- `(a + b) ^ 2^n = a ^ 2^n + b ^ 2^n`. -/
theorem pow_two_pow_add (a b : R) : ∀ n : Nat, (a + b) ^ (2 ^ n) = a ^ (2 ^ n) + b ^ (2 ^ n) := by
  intro n
  induction n with
  | zero => simp [Nat.pow_zero, npow_one]
  | succ n ih =>
      have h2 : 2 ^ (n + 1) = 2 ^ n + 2 ^ n := by rw [Nat.pow_succ, Nat.mul_two]
      rw [h2, pow_add, pow_add, pow_add, ih, sq_add]

/-- The trace `Σ_{i<n} x ^ 2^i`. -/
def frobSum (x : R) : Nat → R
  | 0 => 0
  | n + 1 => frobSum x n + x ^ (2 ^ n)

/-- Telescoping: `Σ_{i<n} ℘(x) ^ 2^i = x ^ 2^n + x`. -/
theorem frobSum_asp (x : R) : ∀ n : Nat, frobSum (asp x) n = x ^ (2 ^ n) + x := by
  intro n
  induction n with
  | zero => simp [frobSum, Nat.pow_zero, npow_one, add_self]
  | succ n ih =>
      have h2 : 2 ^ (n + 1) = 2 ^ n + 2 ^ n := by rw [Nat.pow_succ, Nat.mul_two]
      have hfr : (asp x) ^ (2 ^ n) = (x * x + x) ^ (2 ^ n) := rfl
      rw [frobSum, ih, hfr, pow_two_pow_add, pow_two_pow_mul, h2, pow_add]
      rw [add_comm (x ^ 2 ^ n + x) (x ^ 2 ^ n * x ^ 2 ^ n + x ^ 2 ^ n),
        add_assoc (x ^ 2 ^ n * x ^ 2 ^ n) (x ^ 2 ^ n) (x ^ 2 ^ n + x),
        ← add_assoc (x ^ 2 ^ n) (x ^ 2 ^ n) x, add_self, zero_add]

/-- A field-like characteristic-2 ring on which Frobenius `x ↦ x ^ 2^deg`
is the identity (e.g. F_{2^deg}). -/
class F2FrobField (R : Type u) extends CharTwoRing R where
  deg : Nat
  frob_fix (x : R) : x ^ (2 ^ deg) = x

variable {F : Type u} [F2FrobField F]

/-- **Trace identity**: `Σ_{i<deg} ℘(x) ^ 2^i = 0` over F_{2^deg}: the image of
the Artin–Schreier operator is contained in the trace-zero subspace, the
solubility condition for the AS equations of the GHS tower. -/
theorem trace_asp_zero (x : F) : frobSum (asp x) (F2FrobField.deg F) = 0 := by
  rw [frobSum_asp, F2FrobField.frob_fix, add_self]

end

/-! ### The curve involution of E : y² + xy = x³ + β -/

section Curve
variable {R : Type u} [CharTwoRing R]

/-- Points on the short characteristic-2 curve `y² + xy = x³ + β`. -/
def OnCurve (β x y : R) : Prop := y * y + x * y = x * x * x + β

/-- The map `(x, y) ↦ (x, y + x)` preserves the curve: in the GHS Weil
restriction this involution generates the two irreducible components
(the twist ε ∈ {0,1}) and descends to the hyperelliptic involution of the
descent curve. -/
theorem tau_on (β x y : R) (h : OnCurve β x y) : OnCurve β x (y + x) := by
  unfold OnCurve at h ⊢
  have h1 : (y + x) * (y + x) = y * y + x * x := sq_add y x
  have h2 : x * (y + x) = x * y + x * x := left_distrib x y x
  rw [h1, h2, add_assoc (y * y) (x * x) (x * y + x * x),
    ← add_assoc (x * x) (x * y) (x * x), add_comm (x * x) (x * y),
    add_assoc (x * y) (x * x) (x * x), add_self, add_zero']
  exact h

/-- The involution is an involution: `τ(τ(x, y)) = (x, y)`. -/
theorem tau_tau (x y : R) : ((y + x) + x) = y := by
  rw [add_assoc, add_self, add_zero']

end Curve
