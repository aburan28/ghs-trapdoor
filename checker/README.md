# GHS checker

A dependency-free Python research tool for **classical characteristic-two GHS screening** and **bounded, explicit rational 3-isogeny searches**.

The tool supports ordinary binary elliptic curves

\[
y^2+xy=x^3+ax^2+b,\qquad b\ne0.
\]

It computes exact construction parameters. It does not issue a universal safe/vulnerable verdict, implement arbitrary covers or generalized GHS, construct the descended curve/map, solve the descended DLP, or estimate the complete attack cost. This distinction is part of the JSON report, not just this documentation.

## Run

Python 3.10 or newer is sufficient. No installation, network access, SageMath, Magma, or OpenSSL is needed to run the checker.

```bash
python -m ghs_checker list-presets
python -m ghs_checker check --preset B-233
python -m ghs_checker check --preset teske-secret --subfield-degree 23 --json
python -m ghs_checker field --degree 233 --json
```

Run these commands from the project directory. Optionally install with `python -m pip install .` to obtain the `ghs-checker` command.

The twelve presets include all ten historical NIST B-/K-curves at degrees 163, 233, 283, 409, and 571, plus the explicit secret/public pair from Appendix A of Teske's 2003 preprint. SEC names such as `sect163r2` are accepted as aliases.

## Check a custom curve

```bash
python -m ghs_checker check \
  --degree 7 --modulus 'x^7+x+1' --a 1 --b 0x18 \
  --subgroup-order 19 --json
```

Or use a JSON file:

```json
{
  "name": "my curve",
  "degree": 7,
  "modulus": "x^7+x+1",
  "a": "0x1",
  "b": "0x18",
  "subgroup_order": "19"
}
```

```bash
python -m ghs_checker check --input examples/toy-isogeny.json --json
```

Field elements use a **polynomial basis**: bit `i` is the coefficient of `x^i`. Decimal, `0x` hexadecimal, `0b` binary, and explicit polynomials are accepted. Normal-basis encodings must first be converted; identical hexadecimal strings in two bases do not identify the same element. Values outside the field representation are rejected rather than silently reduced. The defining polynomial is checked for irreducibility.

Optional input fields are `curve_order`, `cofactor`, `source`, `source_url`, and `base_point: {"x": ..., "y": ...}`. A supplied base point is checked against the curve equation. A supplied subgroup order is **claimed input**, not a point-counting or primality certificate. If a curve order is also supplied, divisibility is checked. Reports explicitly leave subgroup-order and descent-image verification unset.

## Explore isogenies

```bash
python -m ghs_checker scan --input examples/toy-isogeny.json \
  --max-curves 8 --max-depth 2 --json
```

The included toy example has 114 points and a subgroup of order 19. The search finds an explicit 3-isogeny from `b=0x18` to `b=0x1e`, lowering the classical parameters from `m=7, genus=63` to `m=4, genus=7`. It also finds another neighboring candidate. This validates **isogeny transport and screening**, not an end-to-end GHS DLP attack or a speedup on this tiny group.

The breadth-first search:

- Enumerates every rational cyclic degree-3 kernel at each expanded node.
- Includes Frobenius-stable kernels whose individual nonzero points are not rational over the input field.
- Records the kernel x-coordinates needed to reconstruct and evaluate the connecting path.
- Checks every proper subfield by default; `--subfield-degree` may be repeated to select subfields.
- Reports curve/depth limits and the fact that other isogeny degrees were not searched.

Defaults are eight curves, depth three, and a screening threshold of genus 32. `--max-genus` changes only this threshold. A curve below it is a candidate for further work, not a verified vulnerability. Runtime grows with the field degree, graph size, and factorization work. Searches never traverse an unbounded graph automatically.

No candidate in a bounded degree-3 search means exactly that. Other prime degrees, generalized descents, alternative covers, and other representatives may still matter. This scanner requires a supplied subgroup order to be coprime to 3. Otherwise injectivity needs point-specific checks, which the scan does not perform.

## Interpret the output

For each descent to `GF(2^s)` the checker returns the extension degree, augmented binary rank `m`, exact genus under the classical regularity condition, relative trace of `b`, and a necessary subgroup-size check.

| Screening label | Meaning |
| --- | --- |
| `small_genus_candidate` | The classical genus is below the configured threshold; a useful map, its kernel, and a practical cost remain to be checked. |
| `above_genus_threshold` | The classical genus exceeds the threshold. This says nothing about arbitrary reductions. |
| `subgroup_size_obstruction` | The Hasse-Weil upper bound is too small to preserve the entire supplied subgroup. |
| `regularity_condition_failed` | The regular classical genus theorem is not asserted for this even-degree/trace case. |

`not_excluded_by_size` is only a necessary condition. It does not prove that a descent map preserves the subgroup. `jacobian_log2_size_leading_term` is `s*g`, the leading size scale `q^g`, **not** a measured group order or attack work factor.

For NIST K-curves the classical genus is one, but the small target over `GF(2)` cannot preserve the large supplied subgroup. For NIST B-curves the direct classical genera are enormous. `field --degree d` separately reproduces the prime-degree Menezes-Qu bounds for `m>1`; these bounds concern that construction, not all possible Weil-descent strategies.

Successful analysis exits with code 0, including reports with no candidate. Invalid or unsupported input exits with code 2. With `--json`, errors are JSON on stderr.

## Python API

```python
from ghs_checker import BinaryField, BinaryCurve, analyze_curve, ThreeIsogeny

field = BinaryField(7, 0x83)  # x^7+x+1
curve = BinaryCurve(field, a=1, b=0x18)
report = analyze_curve(curve, subgroup_order=19)

phi = ThreeIsogeny(curve, kernel_x=0x2)
P = (0x77, 0x11)
assert curve.multiply(19, P) is None
image = phi.apply(P)
assert image is not None
assert phi.target.multiply(19, image) is None
```

## Verification

```bash
python -m unittest discover -v
```

The sixteen tests check the published Teske pair and NIST field bounds, all ten NIST generator orders, both genus branches, the even-degree regularity guard, invariance under Frobenius and an explicit field-basis change, polynomial factorization against exhaustive root enumeration, and explicit isogeny maps against point enumeration, curve orders, and addition. The toy search additionally verifies scalar transport along its reconstructed paths.

Saved JSON examples are in [examples/results](examples/results/). Mathematical definitions and references are in [MATHEMATICS.md](MATHEMATICS.md).

An experimental [genus-2 research extension](research/README.md) audits the exact 2013 GLV/GLS and generic/Kummer parameters and all 45 immediate generic/Kummer Richelot quotient constructions. It also supplies a toy Richelot pair and an explicit genus-4 cover. Its nine additional tests bring the source-tree suite to 25. This extension tests one specified weak-cover class; no new cryptographic vulnerability has been established.

[Alternative descent targets](research/ALTERNATIVE_TARGETS.md) cover direct attacks on Weil restrictions, non-hyperelliptic Jacobians, Pryms, trace-zero varieties, tori, generalized Jacobians and p-adic lifts. A separate [torsion-kernel screener](research/weil_quotient_kernels.py) tests rational polarized quotient kernels at ell=2 or 3 for a cubic elliptic Weil restriction. Nine further tests bring the complete source-tree suite to 34. Positive kernel results still require target and point-map construction.

The [next experiment batch](research/EXPERIMENT_BATCH.md) checks 16,140 actual small cubic-field curve equations, 32 explicit Prym covers, and four torus DLP-transfer examples with smooth elliptic negative controls. It adds eleven tests, bringing the full source-tree suite to **45**. The research modules and saved reports are included in the source tree; the binary-checker wheel does not install them. Positive curve-census results are quotient-kernel candidates, with target reconstruction still pending.
