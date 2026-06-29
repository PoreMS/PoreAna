# v1.0.0

### Performance

Benchmarked on an Apple Silicon macOS machine (2001-frame cylinder trajectory, benzene,
13 lag times) — v1.0.0 local build (`python tests/bench_compare.py`).
Direct comparison against PyPI v0.2.3 is blocked by its CMake/scikit-build dependency,
so microbenchmarks on the core hot paths are provided instead.

| Benchmark | v1.0.0 |
|---|---|
| Density sampling — serial | 0.60 s |
| Density sampling — parallel | 1.33 s |
| Gyration sampling — serial | 0.83 s |
| Diffusion MC sampling — serial | 1.25 s |
| Diffusion MC sampling — parallel | 1.59 s |

Hot-path microbenchmarks (2001 frames, isolated from I/O):

| Operation | Old (Python loops) | New (NumPy) | Speedup |
|---|---|---|---|
| Centre-of-mass per frame | 1.30 s | 0.21 s | **6.2×** |
| Diffusion-bin merge (parallel combine) | 16.3 ms | 4.5 ms | **3.6×** |

Key changes driving the speedup:

* `sample.py` — `_sample_helper` hot loop: position extraction vectorised as `positions[atom_indices] / 10.0 + shift_arr`; COM calculated via `masses_arr @ pos / sum_masses` (NumPy dot product replaces per-atom Python loop); `_masses_arr` pre-computed once in `__init__` as a float64 NumPy array; residue atom index arrays pre-built as `np.array([...])` at init time
* `sample.py` — parallel diffusion-bin merge: replaced O(bin_num × len_window) triple-nested Python loop with `(np.array(a) + np.array(b)).tolist()` per key
* `density.py` — bin volume calculation vectorised: `np.pi * plen * (w[1:]**2 - w[:-1]**2)` replaces scalar loop; weighted-mean integration uses `np.dot` + boolean mask
* `diffusion.py` — MSD normalisation and bin slope calculation fully vectorised with NumPy arrays; Bessel function series uses `np.exp` on a coefficient array; removed unused `math` and `itertools` imports
* `angle.py` / `gyration.py` — density-weighted normalisation uses `np.where(dens != 0, val / dens, 0.0)`; mean computed with `np.mean`; mean line rendered via `np.full_like`
* `adsorption.py` — `np.sum` replaces accumulator loops for reservoir and pore molecule counts
* `mc.py` — `n_proc` parameter no longer shadows the `np` (NumPy) import; duplicate `import numpy as numpy` removed; dead commented-out code removed from `__init__`

### Logic fixes
* `mc.py` — `list_diff_coeff` was incorrectly assigned `list_diff_profile` (per-bin profile data) instead of the actual per-step model coefficients; silently returned wrong values for all MC diffusion runs
* `mc.py` / `sample.py` — function parameter `np=0` shadowed the `numpy` alias, causing `AttributeError: 'int' has no attribute 'array'` in all parallel and post-merge code paths; renamed to `n_proc`
* `diffusion.py` — residual `math.floor` / `math.ceil` calls remained after `import math` was removed; converted to `int(np.floor(...))` / `int(np.ceil(...))`
* `adsorption.py` — reservoir mask `ex_width[:-1]` (150 elements) applied to `ex_bins` (151 elements) caused `IndexError`; mask now uses `ex_width` directly
* `utils.py` / `tables.py` — `file_to_text` and `mc_results` still accessed the old flat pore structure (`pore["box"]`, `pore["diam"]`, `pore["res"]`, `pore["type"]`) after the data model was updated to `pore["box"]["dimensions"]`, `pore["box"]["res"]`, `pore[shape_id]["diam"]`; fixed across all four branches (`gyr_bin`, `diff_bin`, `dens_bin`, `mc`) and in `tables.py`
* `utils.py` — `file_to_text` MC branch used `free_energy[0][i][1:]` (99 elements) in a DataFrame alongside 100-element arrays; `[1:]` removed
* `model.py` — debug `print` calls for `_d0` and `_diff_bin` removed
* `model.py` — duplicate `self._sys_props = {}` assignment removed

### Tests
* Converted `tests/test_simple.py` from unittest to pytest
* Split into `test_unit.py` (fast, no trajectory) and `test_integration.py` (full pipeline)
* Session-scoped `conftest.py` fixture runs all trajectory sampling once per session; pre-computes MC output for `file_to_text` tests
* New coverage: `density.mean()` return values; `adsorption.calculate()` output structure; `density/gyration/angle.bins_plot()` intent "in"/"ex" and normalised x-axis; `diffusion.bins()` output structure; MC output key correctness (`list_diff_coeff` vs `list_diff_profile`); `utils.file_to_text()` for all four output types (dens_bin pore, dens_bin box, diff_bin, gyr_bin, mc); YAML round-trip in `test_utils`
* `bench_compare.py` added: standalone script for sampling speed comparison (`python tests/bench_compare.py`)

### Documentation
* RST source files migrated to MyST Markdown; Sphinx theme updated to furo; API docs via sphinx-autoapi
* Docs source moved from `docsrc/` to `docs/`; previous built HTML preserved at `docs/v_old/`
* Copyright year updated to 2026; DESIGN.md added documenting the red color palette

### CI / tooling
* GitHub Actions: added ruff linting workflow (`lint.yml`)
* GitHub Actions: added pip-audit security scan workflow (`security.yml`)
* CI matrix updated: Python 3.10–3.13; `python_requires` bumped to `>=3.10`

### Administrative
* `setup.py`: version 1.0.0, `python_requires='>=3.10'`, author email updated, chemfiles unpinned, porems dependency bumped to `>=1.0.0`
* README: updated image paths, Python version, PyPI badge, testing and installation instructions


# v0.3.0
* New version due to a change of GitHub organisation.
* You can sample now until a specific frame

# v0.2.3
* New angle routine
* Add option to calculate density in the pore with a constant bin area.
* Add option to specify a upper integration limit when calculating the mean density

# v0.2.2
* Routines can save and load results in hdf5 or obj format
* Add an MC result table to compare diffusion coefficients in self-defined box areas
* Obj. and hdf5 files can be converted into a text file with the most important inputs and results
* Calculate mean density inside Pore
* Restructure sample call for new PoreMS logic - now uses YAML file as input

# v0.2.1
* Most routines are now usable for non-pore simulations
* Add MC routines for calculating diffusion and free energy profiles

# v0.2.0
* Complete restructuring/reprogramming of the code for easier expansion and testing.
