# v1.0.0

### New Features

* **Radial VACF diffusion for cylindrical pores** — `Sample.init_diffusion_vacf()` samples
  the velocity autocorrelation function (VACF) to compute local diffusion coefficients per
  spatial bin. For cylindrical pore systems, `direction="radial_cylindrical"` bins molecules
  radially and decomposes velocities into cylindrical coordinates (radial, tangential, axial).
  Requires a `.trr` trajectory with velocities.
* `diffusion.integrate_bin_diffusion_vacf()` — integrates the sampled VACF using the
  cumulative trapezoid rule to give a time-resolved diffusion coefficient per bin.
* `diffusion.plot_correlation_per_bin()` — visualises the integrated VACF per spatial bin.
* `diffusion.diffusion_per_bin()` — extracts and plots the final diffusion profile.
* `density.density_from_vacf()` / `density_from_vacf_per_residue()` — compute spatial
  density from VACF sampling output.
* `diffusion.mc_profile()` — new `is_legend` parameter to suppress the plot legend.

### Performance

Hot-path microbenchmarks (2001-frame cylinder trajectory, benzene, Apple Silicon macOS,
Python 3.13, isolated from I/O):

| Operation | Old (Python loops) | New (NumPy) | Speedup |
|---|---|---|---|
| Centre-of-mass per frame | 1.30 s | 0.21 s | **6.2x** |
| Diffusion-bin merge (parallel combine) | 16.3 ms | 4.5 ms | **3.6x** |

Key changes driving the speedup:

* `sample.py` — `_sample_helper` hot loop: position extraction vectorised as `positions[atom_indices] / 10.0 + shift_arr`; COM calculated via `masses_arr @ pos / sum_masses` (NumPy dot product replaces per-atom Python loop); `_masses_arr` pre-computed once in `__init__` as a float64 NumPy array; residue atom index arrays pre-built as `np.array([...])` at init time
* `sample.py` — parallel diffusion-bin merge: replaced O(bin_num x len_window) triple-nested Python loop with `(np.array(a) + np.array(b)).tolist()` per key
* `sample.py` / `mc.py` — explicit `mp.get_context("fork")` on non-Windows platforms; Python 3.12 changed the macOS default multiprocessing start method from `fork` to `spawn`, causing each worker to re-import all heavy dependencies on every call — the fix restores ~13x parallel speedup on macOS
* `density.py` — bin volume calculation vectorised: `np.pi * plen * (w[1:]**2 - w[:-1]**2)` replaces scalar loop; weighted-mean integration uses `np.dot` + boolean mask
* `diffusion.py` — MSD normalisation and bin slope calculation fully vectorised with NumPy arrays; Bessel function series uses `np.exp` on a coefficient array
* `angle.py` / `gyration.py` — density-weighted normalisation uses `np.where(dens != 0, val / dens, 0.0)`; mean computed with `np.mean`; mean line rendered via `np.full_like`
* `adsorption.py` — `np.sum` replaces accumulator loops for reservoir and pore molecule counts

### Logic fixes

* `sample.py` — `self.num_res` was never assigned in `__init__`, causing `AttributeError` on all VACF and NumPy batch-processing code paths
* `sample.py` — `in_wall_mask` in `_diffusion_vacf` used the binning-axis coordinate instead of the true radial distance from the pore axis for non-radial pore directions, silently misclassifying wall atoms
* `diffusion.py` — `integrate_bin_diffusion_vacf` divided VACF data by zero for unoccupied bins, producing silent NaN values; unoccupied bins now yield 0.0
* `mc.py` — parallel MC assembly discarded `fluc_diff_bin` / `fluc_df_bin` from all workers except worker 0, so per-bin fluctuation data was silently incomplete after any parallel run
* `diffusion.py` — `mc_profile()` `is_error=True` path accessed `diff_profiles_error_down`/`_up` dicts that were always empty (population code was removed years ago), causing `KeyError`; `is_error` parameter removed and the working `fit.intercept_stderr` error band is now always shown when `infty_profile=True`
* `mc.py` — `list_diff_coeff` was incorrectly assigned `list_diff_profile` (per-bin profile data) instead of the actual per-step model coefficients; silently returned wrong values for all MC diffusion runs
* `mc.py` / `sample.py` — function parameter `np=0` shadowed the `numpy` alias, causing `AttributeError: 'int' has no attribute 'array'` in all parallel and post-merge code paths; renamed to `n_proc`
* `diffusion.py` — residual `math.floor` / `math.ceil` calls remained after `import math` was removed; converted to `int(np.floor(...))` / `int(np.ceil(...))`
* `adsorption.py` — reservoir mask `ex_width[:-1]` (150 elements) applied to `ex_bins` (151 elements) caused `IndexError`; mask now uses `ex_width` directly
* `utils.py` / `tables.py` — `file_to_text` and `mc_results` still accessed the old flat pore structure (`pore["box"]`, `pore["diam"]`, `pore["res"]`, `pore["type"]`) after the data model was updated to `pore["box"]["dimensions"]`, `pore["box"]["res"]`, `pore[shape_id]["diam"]`; fixed across all four output branches
* `utils.py` — `file_to_text` MC branch used `free_energy[0][i][1:]` (99 elements) in a DataFrame alongside 100-element arrays; `[1:]` removed
* `model.py` — debug `print` calls for `_d0` and `_diff_bin` removed; duplicate `self._sys_props = {}` assignment removed

### Removals

* `diffusion.py` — `cui()` (183 lines, zero call sites); radial MC stubs `mc_fit_radial()` and `mc_profile_radial()` preserved as comments for future implementation
* `utils.py` — `column()` and `num_dens_to_mass_dens()` (zero call sites)
* `tables.py` — `mc_statistics()` and `mc_lag_time()` (zero call sites)

### Tests

* Converted `tests/test_simple.py` from unittest to pytest; split into `test_unit.py` (fast, no trajectory) and `test_integration.py` (full pipeline)
* Session-scoped `conftest.py` fixture runs all trajectory sampling once per session; pre-computes MC output for `file_to_text` tests
* New unit tests: `test_vacf_init` (checks `num_res` assignment and `_diffusion_vacf_data` structure); `test_vacf_zero_density` (zero-density bins must yield 0.0, not NaN); `test_diffusion_mc_parallel_keys` (parallel MC must assemble `fluc_diff_bin`/`fluc_df_bin` from all workers)
* New coverage: `density.mean()` return values; `adsorption.calculate()` output structure; `density/gyration/angle.bins_plot()` intent "in"/"ex" and normalised x-axis; `diffusion.bins()` output structure; MC output key correctness (`list_diff_coeff` vs `list_diff_profile`); `utils.file_to_text()` for all four output types; YAML round-trip
* `bench_compare.py` added: standalone script for sampling speed comparison (`python tests/bench_compare.py`)

### Documentation

* Tutorial docs migrated from RST to MyST Markdown; `diffusion_vacf.rst` converted to `diffusion_vacf.md`
* Docs source moved from `docsrc/` to `docs/`; Sphinx theme updated to furo; API docs via sphinx-autoapi
* Copyright year updated to 2026; DESIGN.md added documenting the red color palette

### CI / tooling

* GitHub Actions: added ruff linting workflow (`lint.yml`) and pip-audit security scan (`security.yml`)
* CI matrix updated: Python 3.12-3.13; `python_requires` bumped to `>=3.12`
* Migrated from `setup.py` + `MANIFEST.in` to a single `pyproject.toml` (PEP 517/621); `[tool.ruff]` config added
* CI migrated from `pip` to `uv` (`astral-sh/setup-uv@v5`); `requirements.txt` removed (deps resolved via `pyproject.toml`)
* `codeql.yml`: updated actions from `@v2` to `@v3`; removed irrelevant `javascript` language scan

### Administrative

* `pyproject.toml`: version 1.0.0, `requires-python = ">=3.12"`, author email updated, chemfiles unpinned, porems dependency bumped to `>=1.0.0`; `[project.optional-dependencies]` dev group added (`pytest`, `pytest-cov`)
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
