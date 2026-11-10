"""Sampling speed benchmark: this working tree against a PyPI release.

Builds two throwaway environments with uv, one holding the released poreana
and one holding this tree, runs the same sampling cases in both and prints a
side by side table.

The benchmark bodies are handed to each interpreter with ``python -c``. That
matters for the parallel cases: the worker pools start with *forkserver* or
*spawn*, which re-import the main module, but multiprocessing does not
re-import a ``-c`` main module, so the script is import safe without an
``if __name__ == "__main__":`` guard inside the benchmark body itself.

Usage
-----
    uv run tests/bench_compare.py               # against the latest release
    uv run tests/bench_compare.py --pypi 0.2.2  # against a specific release
    uv run tests/bench_compare.py --quick       # two cases instead of five
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

# Resolve tests/data at import time so the absolute path can be embedded in the
# benchmark script, which runs with a different working directory.
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_BENCH_CASES_FULL = [
    ("density_serial", "density_serial"),
    ("density_parallel", "density_parallel"),
    ("gyration_serial", "gyration_serial"),
    ("diffusion_mc_s", "diffusion_mc_serial"),
    ("diffusion_mc_p", "diffusion_mc_parallel"),
]

_BENCH_CASES_QUICK = [
    ("density_serial", "density_serial"),
    ("gyration_serial", "gyration_serial"),
]


def _make_bench_script(data_dir, out_dir, cases):
    """Return a ``python -c`` body that runs cases and prints JSON results."""
    lines = [
        "import json, os, sys, time, warnings",
        "warnings.filterwarnings('ignore')",
        "import porems as pms",
        "import poreana as pa",
        f"data_dir = {data_dir!r}",
        f"out_dir = {out_dir!r}",
        "mol = pms.Molecule(inp=os.path.join(data_dir, 'benzene.gro'))",
        "traj = os.path.join(data_dir, 'traj_cylinder.xtc')",
        "sys_pore = os.path.join(data_dir, 'pore_system_cylinder_new.yml')",
        "results = {}",
    ]

    # Each body re-creates the Sample, so a case cannot inherit state from the
    # one before it. Every routine writes into out_dir rather than /dev/null,
    # which is not portable and hides a failing save.
    benchmarks = {
        "density_serial": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_density(os.path.join(out_dir, 'dens_s.obj'))",
            "s.sample(is_parallel=False)",
        ],
        "density_parallel": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_density(os.path.join(out_dir, 'dens_p.obj'))",
            "s.sample(is_parallel=True)",
        ],
        "gyration_serial": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_gyration(os.path.join(out_dir, 'gyr_s.obj'))",
            "s.sample(is_parallel=False)",
        ],
        "diffusion_mc_serial": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_diffusion_mc(os.path.join(out_dir, 'mc_s.obj'),"
            " len_step=[1, 2, 5])",
            "s.sample(is_parallel=False)",
        ],
        "diffusion_mc_parallel": [
            "s = pa.Sample(sys_pore, traj, mol)",
            "s.init_diffusion_mc(os.path.join(out_dir, 'mc_p.obj'),"
            " len_step=[1, 2, 5])",
            "s.sample(is_parallel=True)",
        ],
    }

    for label, key in cases:
        lines.append(f"# --- {label} ---")
        lines.append("t0 = time.perf_counter()")
        lines.append("try:")
        for body_line in benchmarks[key]:
            lines.append(f"    {body_line}")
        lines.append(f"    results[{label!r}] = round(time.perf_counter() - t0, 3)")
        lines.append("except Exception as exc:")
        lines.append(f"    results[{label!r}] = 'ERROR: ' + str(exc)")

    lines.append("print(json.dumps(results))")
    return "\n".join(lines)


def _venv_python(venv_dir):
    """Return the interpreter path inside a venv, on POSIX or Windows."""
    if os.name == "nt":
        return os.path.join(venv_dir, "Scripts", "python.exe")
    return os.path.join(venv_dir, "bin", "python")


def _create_env(tmpdir, name, spec):
    """Create a uv venv and install spec into it. Returns its interpreter."""
    venv_dir = os.path.join(tmpdir, name)
    print(f"  creating {venv_dir}", flush=True)
    subprocess.run(["uv", "venv", "--quiet", venv_dir], check=True)

    python_bin = _venv_python(venv_dir)
    print(f"  installing {spec}", flush=True)
    subprocess.run(
        ["uv", "pip", "install", "--quiet", "--python", python_bin, spec],
        check=True,
    )
    return python_bin


def _installed_version(python_bin, cwd):
    """Ask the interpreter which poreana release it ended up with.

    Read from the distribution metadata rather than ``poreana.__version__``,
    which releases before 0.3.0 do not define. ``cwd`` has to be somewhere
    other than the project root: from there the local ``poreana/`` source
    directory and its ``poreana.egg-info`` shadow whatever the venv holds, and
    the probe reports this tree's version for both environments.
    """
    out = subprocess.run(
        [python_bin, "-c",
         "from importlib.metadata import version; print(version('poreana'))"],
        capture_output=True, text=True, cwd=cwd,
    )
    return out.stdout.strip() or "?"


def run_benchmarks(python_bin, script, cwd=None):
    """Run the benchmark body and return the parsed results."""
    result = subprocess.run(
        [python_bin, "-c", script],
        capture_output=True, text=True, cwd=cwd,
    )
    if result.returncode != 0:
        print(f"    stderr: {result.stderr.strip()[:600]}", file=sys.stderr)
        raise RuntimeError(f"benchmark subprocess failed (rc={result.returncode})")
    # The sampling routines print progress, so take the last JSON line.
    for line in reversed(result.stdout.strip().splitlines()):
        if line.startswith("{"):
            return json.loads(line)
    raise ValueError(f"no JSON in output:\n{result.stdout}")


def _fmt(val):
    return f"{val:.3f}s" if isinstance(val, float) else str(val)


def print_table(cases, local_res, pypi_res, pypi_version):
    """Print the side by side table and the geometric mean speedup."""
    labels = [label for label, _ in cases]
    col_w = max(len(label) for label in labels) + 2
    head_local = "this tree"
    head_pypi = f"PyPI {pypi_version}"

    header = (f"{'Benchmark':<{col_w}}  {head_local:>18}  "
              f"{head_pypi:>14}  {'speedup':>8}")
    sep = "-" * len(header)

    print()
    print(sep)
    print(header)
    print(sep)

    for label in labels:
        t_local = local_res.get(label)
        t_pypi = pypi_res.get(label)
        if isinstance(t_local, float) and isinstance(t_pypi, float) and t_local > 0:
            speedup = f"{t_pypi / t_local:.2f}x"
        else:
            speedup = "n/a"
        print(f"{label:<{col_w}}  {_fmt(t_local):>18}  "
              f"{_fmt(t_pypi):>14}  {speedup:>8}")

    print(sep)

    ratios = [
        pypi_res[label] / local_res[label]
        for label in labels
        if isinstance(local_res.get(label), float)
        and isinstance(pypi_res.get(label), float)
        and local_res[label] > 0
    ]
    if ratios:
        gmean = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
        print(f"  geometric mean speedup: {gmean:.2f}x")
    else:
        print("  no case produced a time in both environments")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Compare sampling speed against a released poreana.")
    parser.add_argument("--pypi", metavar="VERSION", default="0.2.3",
                        help="release to compare against (default: 0.2.3)")
    parser.add_argument("--quick", action="store_true",
                        help="run two cases instead of five")
    args = parser.parse_args()

    cases = _BENCH_CASES_QUICK if args.quick else _BENCH_CASES_FULL

    with tempfile.TemporaryDirectory(prefix="poreana_bench_") as tmpdir:
        out_dir = os.path.join(tmpdir, "out")
        os.makedirs(out_dir, exist_ok=True)
        script = _make_bench_script(_DATA_DIR, out_dir, cases)

        print(f"\n[1/2] released poreana {args.pypi}", flush=True)
        pypi_python = _create_env(tmpdir, "pypi_venv", f"poreana=={args.pypi}")
        detected = _installed_version(pypi_python, tmpdir)
        print("  running benchmarks", flush=True)
        pypi_res = run_benchmarks(pypi_python, script, cwd=tmpdir)

        print("\n[2/2] this working tree", flush=True)
        local_python = _create_env(tmpdir, "local_venv", _PROJECT_ROOT)
        print("  running benchmarks", flush=True)
        local_res = run_benchmarks(local_python, script, cwd=tmpdir)

        print_table(cases, local_res, pypi_res, detected)


if __name__ == "__main__":
    main()
