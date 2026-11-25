"""Reproduce the eight finite research tracks and archive their stage costs.

Default: python -m research.run_eight_experiments
Use --sage-command 'docker ... exec CONTAINER sage -python' for an existing
Sage container sharing this checkout. Commands are argv arrays, never a shell.
"""

import argparse
import json
import platform
import shlex
import subprocess
import sys
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
from time import perf_counter, time

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "research/results"


def commands(sage):
    python = [sys.executable, "-m"]
    native = [
        ("wild-classes", "research.wild_artin_schreier"),
        ("wild-classical-comparison", "research.wild_classical_comparison"),
        ("dimension-four-kernels", "research.dimension_four_kernels"),
        ("obstruction-filter", "research.obstruction_filter_probe"),
        ("cost-baseline", "research.eight_cost_baseline"),
    ]
    result = [
        (name, python + [module, "--output", f"research/results/eight-{name}.json"])
        for name, module in native
    ]
    result.append(
        (
            "trace-zero-inputs",
            python
            + [
                "research.quotient_reconstruction",
                "--export-divisor-inputs",
                "--kernel-index",
                "810",
                "--output",
                "research/results/eight-trace-zero-inputs.json",
            ],
        )
    )
    result.append(
        (
            "isogeny-search",
            sage
            + [
                "research/eight_geometry_sage.py",
                "--experiment",
                "1",
                "--output",
                "research/results/eight-isogeny-search.json",
            ],
        )
    )
    result.append(
        (
            "cross-characteristic-targets",
            python
            + [
                "research.cross_characteristic_targets",
                "--output",
                "research/results/eight-cross-characteristic-targets.json",
            ],
        )
    )
    for number, name in (
        (3, "lift-signs"),
        (6, "dimension-four-control"),
        (7, "prym"),
        (8, "wild-cover"),
    ):
        result.append(
            (
                name,
                sage
                + [
                    "research/eight_geometry_sage.py",
                    "--experiment",
                    str(number),
                    "--output",
                    f"research/results/eight-{name}.json",
                ],
            )
        )
    for kind in ("prym", "wild"):
        name = kind + "-subgroup"
        result.append(
            (
                name,
                sage
                + [
                    "research/cover_subgroups_sage.py",
                    "--kind",
                    kind,
                    "--output",
                    f"research/results/eight-{name}.json",
                ],
            )
        )
    result.append(
        (
            "trace-zero-transport",
            sage
            + [
                "research/trace_zero_transport_sage.py",
                "--output",
                "research/results/eight-trace-zero-transport.json",
            ],
        )
    )
    result.append(
        (
            "divisor-inversion-det3",
            sage
            + [
                "research/divisor_inversion_system_sage.py",
                "--budget",
                "300",
                "--output",
                "research/results/eight-divisor-inversion-det3.json",
            ],
        )
    )
    result.append(
        (
            "divisor-inversion-circuit",
            sage
            + [
                "research/divisor_inversion_circuit_sage.py",
                "--budget",
                "120",
                "--output",
                "research/results/eight-divisor-inversion-circuit.json",
            ],
        )
    )
    result.append(
        (
            "independent-checks",
            sage
            + [
                "research/eight_independent_checks_sage.py",
                "--output",
                "research/results/eight-independent-checks.json",
            ],
        )
    )
    result.append(("python-tests", python + ["unittest", "discover", "-v"]))
    return result


def manifest(stages=None):
    circuit = json.loads((RESULTS / "eight-divisor-inversion-circuit.json").read_text())
    if len(circuit["controls"]) != 4:
        raise RuntimeError("The four inversion controls have not completed")
    mapping = {
        "1": ["isogeny-search", "cross-characteristic-targets"],
        "2": ["trace-zero-inputs", "trace-zero-transport"],
        "3": ["lift-signs"],
        "4": [
            "divisor-inversion",
            "divisor-inversion-revised",
            "divisor-inversion-det3",
            "divisor-inversion-circuit",
        ],
        "5": ["obstruction-filter"],
        "6": ["dimension-four-kernels", "dimension-four-control"],
        "7": ["prym", "prym-subgroup"],
        "8": [
            "wild-classes",
            "wild-classical-comparison",
            "wild-cover",
            "wild-subgroup",
        ],
    }
    for names in mapping.values():
        for name in names:
            path = RESULTS / f"eight-{name}.json"
            if not path.is_file():
                raise RuntimeError(f"Missing experiment artifact: {path}")
    artifacts = {}
    for p in sorted(RESULTS.rglob("*")):
        if (
            p.is_file()
            and p.relative_to(RESULTS).parts[0].startswith("eight-")
            and p.name != "eight-manifest.json"
        ):
            artifacts[str(p.relative_to(ROOT))] = {
                "sha256": sha256(p.read_bytes()).hexdigest(),
                "bytes": p.stat().st_size,
            }
    sources = {}
    source_names = [
        "quotient_reconstruction",
        "quotient_divisor_verification",
        "wild_artin_schreier",
        "wild_classical_comparison",
        "dimension_four_kernels",
        "obstruction_filter_probe",
        "eight_geometry_sage",
        "divisor_inversion_system_sage",
        "divisor_inversion_circuit_sage",
        "trace_zero_transport_sage",
        "cover_subgroups_sage",
        "eight_independent_checks_sage",
        "eight_cost_baseline",
        "cross_characteristic_targets",
        "run_eight_experiments",
    ]
    for name in source_names:
        p = ROOT / "research" / (name + ".py")
        sources[str(p.relative_to(ROOT))] = sha256(p.read_bytes()).hexdigest()
    result = {
        "schema_version": 1,
        "base_commit": "85568148d1a04c3c68b0b1c79c2273727cf4bdc8",
        "bounds": "Exhaustive only for the finite matrices in EIGHT_EXPERIMENTS.md; solver timeouts are not impossibility proofs",
        "experiments": mapping,
        "artifacts": artifacts,
        "source_sha256": sources,
        "python_version": platform.python_version(),
        "python_flint_version": version("python-flint"),
        "sage_version": "10.6",
        "Sage_container_digest": "sha256:19995db6194f4a4bab18ce9a88556fd15b9ed5e916b4504fefe618a7796ddbdb",
        "open_goals": [
            "scalable divisor inversion",
            "global normalized lift-sign certificates",
            "normalized chi18 filter",
            "general dimension-four Schottky filter and reconstruction",
            "matching arbitrary restriction factors to Pryms",
            "general ordinary wild-cover subgroup maps",
            "relation generation and linear algebra cost comparison",
        ],
        "practical_speedup_established": False,
    }
    if stages is not None:
        result["stage_costs"] = stages
    (RESULTS / "eight-manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sage-command", default="sage -python")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--manifest-only", action="store_true")
    args = parser.parse_args()
    if args.manifest_only:
        manifest()
        return
    work = commands(shlex.split(args.sage_command))
    if args.list:
        print(json.dumps(work, indent=2))
        return
    logs = RESULTS / "eight-run-logs"
    logs.mkdir(exist_ok=True)
    stages = []
    for name, argv in work:
        start, wall = perf_counter(), time()
        with (logs / (name + ".log")).open("w") as handle:
            result = subprocess.run(
                argv, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, check=False
            )
        stages.append(
            {
                "stage": name,
                "argv": argv,
                "start_unix_seconds": wall,
                "wall_seconds": perf_counter() - start,
                "returncode": result.returncode,
            }
        )
        (RESULTS / "eight-stage-costs.json").write_text(
            json.dumps(stages, indent=2) + "\n"
        )
        if result.returncode:
            raise RuntimeError(f"Stage {name} failed; see {logs / (name + '.log')}")
    manifest(stages)


if __name__ == "__main__":
    main()
