"""
Sampling speed benchmark: local optimized build vs PyPI poreana.

Creates a temporary virtualenv with the latest PyPI release, runs identical
sampling benchmarks in both environments, then prints a side-by-side
comparison table.

Usage
-----
    python tests/bench_compare.py              # compare local vs PyPI 0.2.3
    python tests/bench_compare.py --pypi 0.2.3 # pin a specific PyPI version
    python tests/bench_compare.py --quick      # fewer / smaller cases
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import venv

# ---------------------------------------------------------------------------
# Benchmark suite
#
# The script is sent verbatim to both interpreters via -c.
# All trajectory/system files are referenced by absolute path so the script
# works from any cwd.
# ---------------------------------------------------------------------------

# Resolve the tests/data directory relative to this file at import time so the
# path gets embedded correctly into the benchmark script string.
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

_BENCH_CASES_FULL = [
    ("density_serial",   "density_serial"),
    ("density_parallel", "density_parallel"),
    ("gyration_serial",  "gyration_serial"),
    ("diffusion_mc_s",   "diffusion_mc_serial"),
    ("diffusion_mc_p",   "diffusion_mc_parallel"),
]

_BENCH_CASES_QUICK = [
    ("density_serial",  "density_serial"),
    ("gyration_serial", "gyration_serial"),
]


def _make_bench_script(data_dir, cases):
    """Return a Python -c script that runs cases and prints JSON results."""
    lines = [
        "import json, time, sys, warnings, os",
        "warnings.filterwarnings('ignore')",
        "import porems as pms",
        "import poreana as pa",
        f"data_dir = {data_dir!r}",
        "print('poreana version: ' + getattr(pa, '__version__', '?'), file=sys.stderr)",
        # Warm-up: load modules so lazy imports don't count toward timing.
        "try:",
        "    mol = pms.Molecule(inp=os.path.join(data_dir, 'benzene.gro'))",
        "    s = pa.Sample(os.path.join(data_dir, 'pore_system_cylinder_new.yml'),",
        "                  os.path.join(data_dir, 'traj_cylinder.xtc'), mol)",
        "except Exception as _e:",
        "    print('warm-up failed:', _e, file=sys.stderr)",
        "results = {}",
    ]

    setup = [
        "mol = pms.Molecule(inp=os.path.join(data_dir, 'benzene.gro'))",
        "traj = os.path.join(data_dir, 'traj_cylinder.xtc')",
        "sys_pore = os.path.join(data_dir, 'pore_system_cylinder_new.yml')",
    ]

    benchmarks = {
        "density_serial": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_density('/dev/null')",
            "s.sample(is_parallel=False)",
        ],
        "density_parallel": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_density('/dev/null')",
            "s.sample(is_parallel=True)",
        ],
        "gyration_serial": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_gyration('/dev/null')",
            "s.sample(is_parallel=False)",
        ],
        "diffusion_mc_serial": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_diffusion_mc('/tmp/_bench_mc.obj', len_step=[1, 2, 5])",
            "s.sample(is_parallel=False)",
        ],
        "diffusion_mc_parallel": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_diffusion_mc('/tmp/_bench_mc_p.obj', len_step=[1, 2, 5])",
            "s.sample(is_parallel=True)",
        ],
    }

    case_keys = {label: key for label, key in cases}

    for label, key in cases:
        bench_lines = benchmarks.get(key, [])
        lines.append(f"# --- {label} ---")
        for ln in setup:
            lines.append(ln)
        lines.append(f"t0 = time.perf_counter()")
        lines.append(f"try:")
        for ln in bench_lines:
            lines.append(f"    {ln}")
        lines.append(f"    results[{label!r}] = round(time.perf_counter() - t0, 3)")
        lines.append(f"except Exception as _e:")
        lines.append(f"    results[{label!r}] = f'ERROR: {{_e}}'")

    lines.append("print(json.dumps(results))")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Venv helpers
# ---------------------------------------------------------------------------

def _create_pypi_venv(version, tmpdir):
    """Create a venv and pip-install poreana==version from PyPI."""
    venv_dir = os.path.join(tmpdir, "pypi_venv")
    print(f"  Creating venv at {venv_dir} ...", flush=True)
    venv.create(venv_dir, with_pip=True, clear=True)

    pip = os.path.join(venv_dir, "bin", "pip")
    spec = f"poreana=={version}"
    print(f"  Installing {spec} from PyPI ...", flush=True)
    subprocess.run(
        [pip, "install", "--quiet", "--no-cache-dir", spec, "scipy", "chemfiles", "porems"],
        check=True,
        cwd=tmpdir,
    )
    return os.path.join(venv_dir, "bin", "python")


def _create_local_venv(tmpdir):
    """Create a venv and install the local poreana source tree."""
    venv_dir = os.path.join(tmpdir, "local_venv")
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print(f"  Creating venv at {venv_dir} ...", flush=True)
    venv.create(venv_dir, with_pip=True, clear=True)

    pip = os.path.join(venv_dir, "bin", "pip")
    print(f"  Installing local build from {project_root} ...", flush=True)
    subprocess.run(
        [pip, "install", "--quiet", project_root],
        check=True,
    )
    return os.path.join(venv_dir, "bin", "python")


# ---------------------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------------------

def run_benchmarks(python_bin, script, cwd=None):
    """Run bench script; return dict of results."""
    result = subprocess.run(
        [python_bin, "-c", script],
        capture_output=True,
        text=True,
        cwd=cwd,
    )
    if result.returncode != 0:
        print(f"    stderr: {result.stderr.strip()[:600]}", file=sys.stderr)
        raise RuntimeError(f"Benchmark subprocess failed (rc={result.returncode})")
    for line in reversed(result.stdout.strip().splitlines()):
        if line.startswith("{"):
            return json.loads(line)
    raise ValueError(f"No JSON found in output:\n{result.stdout}")


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _fmt(val):
    if isinstance(val, float):
        return f"{val:.3f}s"
    return str(val)


def print_table(cases, local_res, pypi_res, pypi_version):
    import math
    labels = [c[0] for c in cases]
    col_w = max(len(lb) for lb in labels) + 2
    ver_local = "local (optimized)"
    ver_pypi  = f"PyPI {pypi_version}"

    header = f"{'Benchmark':<{col_w}}  {ver_local:>18}  {ver_pypi:>14}  {'speedup':>8}"
    sep = "-" * len(header)

    print()
    print(sep)
    print(header)
    print(sep)

    for label in labels:
        t_local = local_res.get(label)
        t_pypi  = pypi_res.get(label)

        if isinstance(t_local, float) and isinstance(t_pypi, float) and t_pypi > 0:
            speedup = f"{t_pypi / t_local:.2f}×"
            marker  = "  ✓" if t_pypi > t_local else "  ="
        else:
            speedup = "n/a"
            marker  = ""

        print(f"{label:<{col_w}}  {_fmt(t_local):>18}  {_fmt(t_pypi):>14}  {speedup:>8}{marker}")

    print(sep)

    ratios = [
        pypi_res[lb] / local_res[lb]
        for lb in labels
        if isinstance(local_res.get(lb), float)
        and isinstance(pypi_res.get(lb), float)
        and local_res[lb] > 0
    ]
    if ratios:
        gmean = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
        print(f"  Geometric mean speedup: {gmean:.2f}×")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pypi", metavar="VERSION", default="0.2.3",
                        help="PyPI version to compare against (default: 0.2.3)")
    parser.add_argument("--quick", action="store_true",
                        help="Run a smaller subset of benchmarks")
    args = parser.parse_args()

    cases = _BENCH_CASES_QUICK if args.quick else _BENCH_CASES_FULL
    script = _make_bench_script(_DATA_DIR, cases)

    with tempfile.TemporaryDirectory(prefix="poreana_bench_") as tmpdir:
        # ---- PyPI version ------------------------------------------------
        print(f"\n[1/2] Installing PyPI poreana ({args.pypi}) ...", flush=True)
        pypi_python = _create_pypi_venv(args.pypi, tmpdir)

        print(f"\n[1/2] Running benchmarks with PyPI version ...", flush=True)
        pypi_res = run_benchmarks(pypi_python, script, cwd=tmpdir)

        pip_pypi = os.path.join(os.path.dirname(pypi_python), "pip")
        ver_out = subprocess.run(
            [pip_pypi, "show", "poreana"],
            capture_output=True, text=True, cwd=tmpdir,
        )
        detected = args.pypi
        for line in ver_out.stdout.splitlines():
            if line.startswith("Version:"):
                detected = line.split(":", 1)[1].strip()
                break

        # ---- Local version -----------------------------------------------
        print(f"\n[2/2] Installing local build ...", flush=True)
        local_python = _create_local_venv(tmpdir)

        print(f"\n[2/2] Running benchmarks with local build ...", flush=True)
        local_res = run_benchmarks(local_python, script, cwd=tmpdir)

        # ---- Report -------------------------------------------------------
        print_table(cases, local_res, pypi_res, detected)


if __name__ == "__main__":
    main()
