# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `pyproject.toml` with build system and dependencies.
- `CITATION.cff` which drives GitHub's "Cite this repository" button.
- Local diffusion analysis based on the velocity autocorrelation function.
  `Sample.init_diffusion_vacf` samples the VACF along a box axis or, inside a
  cylindrical pore, in radial, axial and tangential direction;
  `diffusion.plot_correlation_per_bin`, `diffusion.diffusion_per_bin` and
  `diffusion.integrate_bin_diffusion_vacf` evaluate the sampled file.

### Changed

- Documentation sources moved from `docsrc/` to `docs/source/`, Bootstrap theme
  to furo, and converted from `.rst` to MyST `.md`
- The sampling loop now evaluates whole frames with numpy masks instead of
  iterating over single atoms.
- **Breaking:** `Sample.sample` renamed its `np` argument to `num_cores` and
  `MC.run` takes `n_proc` instead of `np`.
- **Breaking:** parallel routines start their workers with `forkserver` on POSIX and
  `spawn` on Windows instead of `fork`. Therefore, a script that runs a parallel
  should be guarded with `if __name__ == "__main__":`.
- The test suite moved from unittest to pytest with more and better tests.

### Fixed

- A series of bugs, including reading sampling files, parallel data
  collection, frame reshaping, and numerical edge cases.

### Removed

- `setup.py`, the empty `MANIFEST.in`, `requirements.txt`,
  `docsrc/requirements.txt`, `docsrc/html.sh` and `tests/venv.txt`.
- **Breaking:** hdf5 is no longer supported.
- **Breaking:** `tables.mc_statistics` and `tables.mc_lag_time`, which nothing used.
- **Breaking:** `MC.run(do_radial=...)` with unfinished radial diffusion code

## [0.3.0] - 2025-10-23

### Added

- Sampling can now be stopped at a specific frame.

### Changed

- New version due to a change of GitHub organisation.

## [0.2.3] - 2024-11-08

### Added

- New angle routine.
- Option to calculate the density in the pore with a constant bin area.
- Option to specify an upper integration limit when calculating the mean
  density.

## [0.2.2] - 2022-02-03

### Added

- Routines can save and load results in hdf5 or obj format.
- An MC result table to compare diffusion coefficients in self-defined box
  areas.
- Obj and hdf5 files can be converted into a text file holding the most
  important inputs and results.
- Calculation of the mean density inside a pore.

### Changed

- The sample call was restructured for the new PoreMS logic and now takes a
  YAML file as input.

## [0.2.1] - 2021-06-29

### Added

- MC routines for calculating diffusion and free energy profiles.

### Changed

- Most routines are now usable for non-pore simulations.

## [0.2.0] - 2021-05-05

### Changed

- Complete restructuring and reprogramming of the code for easier expansion
  and testing.

[0.2.3]: https://pypi.org/project/PoreAna/0.2.3/
[0.2.2]: https://pypi.org/project/PoreAna/0.2.2/
[0.2.1]: https://pypi.org/project/PoreAna/0.2.1/
[0.2.0]: https://pypi.org/project/PoreAna/0.2.0/
