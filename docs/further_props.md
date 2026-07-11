# Further Properties

Other properties can be implemented using available analysis
routines as a basis. Following analyses have also been implemented.

## Adsorption

Based on the Density analysis

```python
import porems as pms
import poreana as pa

mol = pms.Molecule(inp="data/benzene.gro")

sample = pa.Sample("data/pore_system_cylinder.obj", "data/traj_cylinder.xtc", mol, [])
sample.init_density("output/dens.obj")
sample.sample()
```


`Finished frame 2001/2001...`


Adsorption values within the pore can be calculated by counting the
instance of atoms inside and outside the pore system resulting into a
surface concentration within the pore and a volumetric concentration in
the bulk reservoirs

```python
import pandas as pd

ads = pa.adsorption.calculate("output/dens.obj")

pd.DataFrame(ads)
```


## Radius of Gyration

The radius of gyration can be calculated by sampling the molecule radius
inside and outside the pore using the gyration routine. The density
routine is also needed for weighting the radius based on the bin
density

```python
import porems as pms
import poreana as pa

mol = pms.Molecule(inp="data/benzene.gro")

sample = pa.Sample("data/pore_system_cylinder.obj", "data/traj_cylinder.xtc", mol, [])
sample.init_density("output/dens.obj")
sample.init_gyration("output/gyr.obj")
sample.sample()
```


`Finished frame 2001/2001...`


The gyration radius can then be calculated and visualized as a function
of radius and distance inside and outside the pore respectively

```python
import matplotlib.pyplot as plt

plt.figure(figsize=(13, 4))

ylim = [0, 0.2]

plt.subplot(121)
pa.gyration.bins_plot("output/gyr.obj", "output/dens.obj", intent="in")
plt.xlim([0, 2])
plt.ylim(ylim)
plt.xlabel("Distance from pore center (nm)")
plt.ylabel(r"Radius of gyration (nm)")

plt.subplot(122)
pa.gyration.bins_plot("output/gyr.obj", "output/dens.obj", intent="ex")
plt.xlim([0, 5])
plt.ylim(ylim)
plt.xlabel("Distance from reservoir end (nm)")
plt.ylabel(r"Radius of gyration (nm)")
```


:::{figure} /pics/gyration_01.svg
:align: center
:width: 100%
:::


## Angle

The angle can be calculated by sampling the angles between a molecule vector,
defined by two atom ids, and the surface normal vector inside and outside the
pore using the angle routine. The density routine is also needed for weighting
the angle based on the bin density

```python
import porems as pms
import poreana as pa

mol = pms.Molecule(inp="data/benzene.gro")

sample = pa.Sample("data/pore_system_cylinder.obj", "data/traj_cylinder.xtc", mol, [])
sample.init_density("output/dens.obj")
sample.init_angle("output/angle.obj", [0, 3])
sample.sample(is_parallel=False)
```


`Finished frame 2001/2001...`


Note that the angle routine cannot be parallelized. Optionally, or if the pore
shape has not been implemented, the PoreMS *Shape* module can be used for
determining the normal vectors

```python
shape = pms.Cylinder({"centroid": centroid, "central": [0, 0, 1], "length": length, "diameter": diameter})
def normal_in(pos): return shape.normal(pos)
def normal_ex(pos): return [0, 0, -1] if pos[2] < centroid[2] else [0, 0, 1]
normals = {"in": normal_in, "ex": normal_ex}

sample.init_angle("output/angle.obj", [0, 3], normals=normals)
```


For more information on shapes, visit the [PoreMS documentation](https://porems.github.io/PoreMS/autoapi/index.html).
The angles can then be calculated and visualized as a function
of distance inside and outside the pore respectively

```python
import matplotlib.pyplot as plt

plt.figure(figsize=(13, 4))

ylim = [0, 0.2]

plt.subplot(121)
pa.angle.bins_plot("output/angle.obj", "output/dens.obj", intent="in")
plt.xlim([0, 2])
plt.ylim(ylim)
plt.xlabel("Distance from pore center (nm)")
plt.ylabel(r"Angle (deg)")

plt.subplot(122)
pa.angle.bins_plot("output/angle.obj", "output/dens.obj", intent="ex")
plt.xlim([0, 5])
plt.ylim(ylim)
plt.xlabel("Distance from reservoir end (nm)")
plt.ylabel(r"Angle (deg)")
```


:::{figure} /pics/angle_01.svg
:align: center
:width: 100%
:::
