# Diffusion Analysis in a Pore (Bin Method)

The bin diffusion analysis needs the sampled object file using the bin
diffusion routine

```python
import porems as pms
import poreana as pa

mol = pms.Molecule(inp="data/benzene.gro")

sample = pa.Sample("data/pore_system_cylinder.obj", "data/traj_cylinder.xtc", mol, [])
sample.init_diffusion_bin("output/diff.obj")
sample.sample()
```

`Finished frame 2001/2001...`


The diffusion can be calculated through the entire pore

```python
pa.diffusion.cui("output/diff.obj")
```

`Diffusion axial:  0.741 10^-9 m^2s^-1`

`Diffusion radial: 648.838 10^-9 m^2 s^-1; Number of zeros: 41; Radius:  0.17`

:::{figure} /pics/diffusion_bins_01.svg
:align: center
:width: 50%
:::


or binwise

```python
pa.diffusion.bins("output/diff.obj")
```


:::{figure} /pics/diffusion_bins_02.svg
:align: center
:width: 50%
:::


By weighting the diffusion profile with the density profile, one can
calculate an axial mean value

```python
pa.diffusion.mean("output/diff.obj", "output/dens.obj")
```


`Mean Diffusion axial: 1.131 10^-9 m^2s^-1`
