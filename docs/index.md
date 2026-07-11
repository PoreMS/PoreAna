:::{figure} /pics/logo_text.svg
:align: center
:width: 40%
:::

# Pore System Analysis Tool

:::{note}
You are viewing the current documentation. The archived documentation for
the previous release is available at [v_old docs](v_old/index.html).
:::

This Python package analyzes GROMACS pore simulations produced with
[PoreSim](https://porems.github.io/PoreSim/).
To see the code or report a bug, please visit the
[GitHub repository](https://github.com/PoreMS/PoreAna).

## Tutorials

- [Sampling](sample.md) — trajectory sampling setup
- [Density](density.md) — radial and axial density profiles
- [Diffusion (Binning)](diffusion_bin.md) — mean-squared displacement via binning
- [Diffusion (MC)](diffusion_mc.md) — Monte Carlo diffusion analysis
- [Diffusion (VACF)](diffusion_vacf.md) — velocity autocorrelation function diffusion
- [Further Properties](further_props.md) — gyration radius, angles, adsorption

## API Reference

Full API documentation is auto-generated from source: [API Reference](autoapi/index.html)

:::{toctree}
:hidden:
:maxdepth: 2

sample
density
diffusion_bin
diffusion_mc
diffusion_vacf
further_props
autoapi/index
:::
