---
hide-toc: true
---

# PoreAna

PoreAna analyzes molecular-dynamics trajectories of nanopore systems and plain
simulation boxes. Trajectories are read with
[chemfiles](https://pypi.org/project/chemfiles/), the pore geometry comes from a
[PoreMS](https://porems.github.io/PoreMS/) system file, and the analysis routines
produce density, radius-of-gyration, angle, diffusion and free-energy profiles.

To see the code or report a bug, please visit the
[GitHub repository](https://github.com/PoreMS/PoreAna).

## Typical workflow

1. [Sample the trajectory](sample.md) once, activating the analyses you need.
2. Evaluate the sampled data:
   - [Density analysis](density.md)
   - [Diffusion, bin method](diffusion_bin.md)
   - [Diffusion, Monte-Carlo method](diffusion_mc.md)
   - [Adsorption, radius of gyration and angle](further_props.md)

The [API reference](api.rst) documents every module in full.

## Installation

```bash
pip install poreana
```

```{toctree}
:hidden:

sample
density
diffusion_bin
diffusion_mc
further_props
```

```{toctree}
:hidden:
:caption: API

api
```
