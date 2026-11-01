---
hide-toc: true
---

# PoreAna

```{include} ../../README.md
:start-after: <!-- start elevator-pitch -->
:end-before: <!-- end elevator-pitch -->
```

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

```{include} ../../README.md
:start-after: <!-- start installation -->
:end-before: <!-- end installation -->
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
