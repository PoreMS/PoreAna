<img src="https://github.com/PoreMS/PoreAna/blob/main/docs/pics/logo_text_sub.svg" width="60%">

--------------------------------------

[![PyPI Version](https://img.shields.io/badge/PyPI-1.0.0-orange)](https://pypi.org/project/poreana/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/PoreMS/PoreAna/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14056630.svg)](https://doi.org/10.5281/zenodo.14056630)
[![Build Status](https://github.com/PoreMS/PoreAna/actions/workflows/workflow.yml/badge.svg)](https://github.com/PoreMS/PoreAna/actions/workflows/workflow.yml)
[![codecov](https://codecov.io/gh/PoreMS/PoreAna/branch/main/graph/badge.svg)](https://codecov.io/gh/PoreMS/PoreAna)

## Documentation

Online documentation is available at [porems.github.io/PoreAna](https://porems.github.io/PoreAna/).

The docs include examples for [density analysis](https://porems.github.io/PoreAna/density.html),
[diffusion](https://porems.github.io/PoreAna/diffusion_mc.html), and an [API reference](https://porems.github.io/PoreAna/autoapi/index.html).


## Dependencies

PoreAna requires Python 3.12+.

Installation requires [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/),
[chemfiles](https://pypi.org/project/chemfiles/), [seaborn](https://seaborn.pydata.org/),
[scipy](https://scipy.org/), [h5py](https://pypi.org/project/h5py/), and
[porems](https://pypi.org/project/porems/) >= 1.0.0.


## Installation

The latest stable release can be installed from PyPI:

    pip install poreana

Or install the development version directly from GitHub:

    pip install git+https://github.com/PoreMS/PoreAna.git#egg=poreana

Or download the repository and install in the top directory via:

    pip install .


## Testing

Install in editable mode with test dependencies:

    pip install -e ".[dev]"

Then run the tests:

    pytest tests/test_unit.py          # fast unit tests
    pytest tests/test_integration.py   # full integration tests (requires data files, slow)


## Development

PoreAna development takes place on GitHub: [www.github.com/PoreMS/PoreAna](https://github.com/PoreMS/PoreAna)

Please submit any reproducible bugs you encounter to the [issue tracker](https://github.com/PoreMS/PoreAna/issues).


## How to Cite PoreAna

When citing PoreAna please use the current **Zenodo DOI** corresponding to the used PoreAna version. (Current DOI is listed in the badges.)
