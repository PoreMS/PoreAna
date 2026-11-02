"""Fast tests that do not run a full sampling pass."""

import numpy as np
import pytest
import yaml

import porems as pms
import poreana as pa


#########
# Utils #
#########
def test_utils():
    pa.utils.mkdirp("output/test")

    # Row to column matrix
    assert pa.utils.column([[1, 1, 1], [2, 2, 2]]) == [[1, 2], [1, 2], [1, 2]]

    # obj round-trip
    file_link = "output/test/test.obj"
    pa.utils.save([1, 1, 1], file_link)
    assert pa.utils.load(file_link) == [1, 1, 1]

    # Unknown extension is rejected on both ends
    pa.utils.save(file_link, "data/data.DOTA")
    assert pa.utils.load("data/data.DOTA") is None
    assert pa.utils.load("output/test/test.txt", file_type="DOTA") is None

    # Unit conversions
    assert round(pa.utils.mumol_m2_to_mols(3, 100), 4) == 180.66
    assert round(pa.utils.mols_to_mumol_m2(180, 100), 4) == 2.989
    assert round(pa.utils.mmol_g_to_mumol_m2(0.072, 512), 2) == 0.14
    assert round(pa.utils.mmol_l_to_mols(30, 1000), 4) == 18.066
    assert round(pa.utils.mols_to_mmol_l(18, 1000), 4) == 29.8904

    pa.utils.toc(pa.utils.tic(), message="Test", is_print=True)
    assert round(pa.utils.toc(pa.utils.tic(), is_print=True)) == 0


@pytest.mark.xfail(strict=True, reason="numpy types are dumped as python/object/apply tags")
def test_utils_yaml_is_plain():
    """A saved .yml must be readable without yaml's unsafe loader."""
    pa.utils.mkdirp("output/test")
    link = "output/test/numpy.yml"
    pa.utils.save({"a": np.float64(1.5), "b": np.arange(3), "c": np.int64(7)}, link)

    with open(link) as f:
        plain = yaml.safe_load(f)
    assert plain == {"a": 1.5, "b": [0, 1, 2], "c": 7}

    loaded = pa.utils.load(link)
    assert type(loaded["a"]) is float
    assert type(loaded["b"]) is list
    assert type(loaded["c"]) is int


############
# Geometry #
############
def test_geometry():
    vec_a = [1, 1, 2]
    vec_b = [0, 3, 2]
    assert round(pa.geom.dot_product(vec_a, vec_b), 4) == 7
    assert round(pa.geom.length(vec_a), 4) == 2.4495
    assert [round(x, 4) for x in pa.geom.vector(vec_a, vec_b)] == [-1, 2, 0]
    assert pa.geom.vector([0, 1], [0, 0, 0]) is None
    assert [round(x, 4) for x in pa.geom.unit(vec_a)] == [0.4082, 0.4082, 0.8165]
    assert [round(x, 4) for x in pa.geom.cross_product(vec_a, vec_b)] == [-4, -2, 3]
    assert round(pa.geom.angle(vec_a, vec_b), 4) == 37.5714


def test_geometry_angle_parallel():
    """Parallel and antiparallel vectors must not raise a math domain error.

    dot(a,b)/(|a||b|) overshoots 1 by one ulp for a large share of vectors -
    about a quarter of random triples - and acos then rejects it.
    """
    vec = [1, 1, 2]
    assert round(pa.geom.angle(vec, vec), 4) == 0.0
    assert round(pa.geom.angle(vec, [-x for x in vec]), 4) == 180.0
    assert round(pa.geom.angle([3, 1, 4], [3, 1, 4]), 4) == 0.0


###########
# Density #
###########
def test_density_mean_integrates_every_bin():
    """density.mean must integrate all radial bins, as its docstring states.

    The loop runs over range(bin_num-1), so the outermost bin never
    contributes. It stays invisible whenever that bin is empty, which it is
    for the bundled trajectories.
    """
    density = {
        "sample": {
            "pore": {"shape_00": {}},
            "inp": {"mass": 1.0},
            "data": {"shape_00": {"in_width": [0.0, 1.0, 2.0, 3.0, 4.0]}},
        },
        "num_dens": {"shape_00": {"in": [4.0, 3.0, 2.0, 1.0]}},
    }

    # (4*1 + 3*3 + 2*5 + 1*7) / (1 + 3 + 5 + 7)
    mean = pa.density.mean(density, is_print=False, int_limit=10)
    assert round(mean["num_dens"]["shape_00"], 6) == 1.875


##############
# Pore shapes #
##############
def test_bin_in_const_A_per_shape(change_to_tests_dir):
    """Every shape must get the same number of constant-area bins.

    bin_num is incremented inside the shape loop, so the second shape is
    binned with bin_num+4 and the third with bin_num+6.
    """
    mol = pms.Molecule(inp="data/benzene.gro")
    sample = pa.Sample("data/pore.yml", "data/traj_cylinder.xtc", mol)

    result = sample._bin_in_const_A(10)
    counts = {pore_id: len(w) for pore_id, w in result["width"].items()}
    assert set(counts) == {"shape_00", "shape_01"}
    assert counts["shape_00"] == counts["shape_01"] == 12
    assert len(result["bins"]) == 13


def test_init_angle_normals_per_shape(change_to_tests_dir):
    """Each shape must get its own surface normal function.

    normal_in is defined inside the shape loop and closes over `shape`, so
    every entry ends up evaluating the last shape's normal.
    """
    mol = pms.Molecule(inp="data/benzene.gro")
    sample = pa.Sample("data/pore.yml", "data/traj_cylinder.xtc", mol)
    sample.init_angle("output/test/angle_two_shapes.obj", [0, 3])

    normals = sample._angle_normals
    pos = [4.0, 4.0, 5.0]
    # shape_00 is centred on [3, 3, 5], shape_01 on [8, 8, 5]
    assert [round(float(x), 4) for x in normals["shape_00"]["in"](pos)] == [1.0, 1.0, 0.0]
    assert [round(float(x), 4) for x in normals["shape_01"]["in"](pos)] == [-4.0, -4.0, 0.0]


################
# MC diffusion #
################
def test_init_diffusion_mc_keeps_len_step(change_to_tests_dir):
    """init_diffusion_mc must not reorder the list it was handed."""
    mol = pms.Molecule(inp="data/benzene.gro")
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)

    len_step = [10, 2, 5]
    sample.init_diffusion_mc("output/test/mc.obj", len_step=len_step)
    assert len_step == [10, 2, 5]


#########
# Model #
#########
def _write_minimal_mc_sample(link):
    """Write the smallest sample file a Model will accept."""
    pa.utils.mkdirp("output/test")
    bin_num = 4
    pa.utils.save({
        "type": "diff_mc",
        "box": {"length": [1.0, 1.0, 4.0]},
        "inp": {"bin_num": bin_num, "num_frame": 10, "len_step": [1],
                "len_frame": 2e-12, "bins": [0.0, 1.0, 2.0, 3.0, 4.0],
                "direction": 2, "is_pbc": True},
        "data": {1: np.zeros((bin_num+2, bin_num+2))},
    }, link)
    return link


def test_model_accepts_float():
    link = _write_minimal_mc_sample("output/test/mc_minimal.obj")
    model = pa.CosineModel(link, 2, 2, d0=1.0)
    # _d0 = 1.0 * 1e9/1e12, bin width 1 nm -> diff_unit = log(1) = 0
    assert np.array_equal(np.round(model._diff_bin, 4), np.array([-6.9078] * model._bin_num))


def test_model_accepts_numpy_float():
    """A numpy float guess must initialize the diffusion profile.

    The guard is `type(d0) == float`, which no numpy scalar satisfies, so
    self._d0 is never assigned at all.
    """
    link = _write_minimal_mc_sample("output/test/mc_minimal.obj")
    model = pa.CosineModel(link, 2, 2, d0=np.float64(1.0))
    assert np.array_equal(np.round(model._diff_bin, 4), np.array([-6.9078] * model._bin_num))


def test_model_accepts_int_d0():
    """An integer guess must give the same profile as the equal float.

    StepModel's own default is the integer 1, so an integer is a natural
    thing to hand CosineModel, but isinstance(1, float) is False.
    """
    link = _write_minimal_mc_sample("output/test/mc_minimal.obj")
    reference = pa.CosineModel(link, 2, 2, d0=1.0)
    model = pa.CosineModel(link, 2, 2, d0=1)
    assert np.array_equal(model._diff_bin, reference._diff_bin)
    assert np.array_equal(model._diff_coeff, reference._diff_coeff)


def test_model_rejects_list_d0_without_coeff():
    """A per-bin guess needs matching coefficients, and must say so.

    The list branch assigns self._coeff to the diffusion coefficients, so an
    omitted coeff used to fail far away in _calc_profile with
    "operands could not be broadcast together with shapes (0,) (4,2)".
    """
    link = _write_minimal_mc_sample("output/test/mc_minimal.obj")
    with pytest.raises(ValueError, match="coeff"):
        pa.CosineModel(link, 2, 2, d0=[1.0, 1.0, 1.0, 1.0])

    # A matching coeff is accepted
    model = pa.CosineModel(link, 2, 2, d0=[1.0, 1.0, 1.0, 1.0], coeff=np.zeros(2))
    assert isinstance(model._d0, list)

    # And an outright unusable type is named rather than causing an
    # AttributeError on self._d0 later
    with pytest.raises(ValueError, match="d0"):
        pa.CosineModel(link, 2, 2, d0="1.0")


########
# VACF #
########
def test_vacf_data_structure(change_to_tests_dir):
    """_bin_pore and _diffusion_vacf_data must agree on the bin count."""
    mol = pms.Molecule(inp="data/benzene.gro")
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)

    # Radial bin edges of a cylindrical pore
    bins = sample._bin_pore(16)["bins"]
    for pore_id, edges in bins.items():
        assert pore_id.startswith("shape")
        assert edges[0] == 0.0
        assert len(edges) == 17

    assert sample.num_res == len(sample._res_list)

    sample._is_diffusion_vacf = True
    sample._diff_vacf_inp = {
        "bin_num": 8,
        "corr_steps": 10,
        "sample_each_residue": False,
        "direction": 2,
        "bins": {"shape_00": np.linspace(0, 1, 9)},
    }
    data = sample._diffusion_vacf_data()
    assert set(data) == {"shape_00"}
    # bins x corr_steps x residues x directions
    assert data["shape_00"]["vacf_data"].shape == (8, 10, 1, 3)
    assert data["shape_00"]["density"].shape == (8, 1)


def test_init_diffusion_vacf_rejects_xtc(change_to_tests_dir):
    """VACF needs velocities, so an .xtc trajectory must be refused."""
    mol = pms.Molecule(inp="data/benzene.gro")
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_diffusion_vacf("output/test/vacf.obj")
    assert not sample._is_diffusion_vacf


def test_init_diffusion_vacf_rejects_bin_conflict(change_to_tests_dir):
    """Bin diffusion and VACF cannot be sampled in the same run."""
    mol = pms.Molecule(inp="data/benzene.gro")
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample._traj = "dummy.trr"
    sample._is_diffusion_bin = True
    sample.init_diffusion_vacf("output/test/vacf.obj", direction=2)
    assert not sample._is_diffusion_vacf


def _vacf_sample(bin_num=4, direction="radial"):
    """A VACF sample on the single-cylinder pore, without a trajectory."""
    mol = pms.Molecule(inp="data/benzene.gro")
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample._traj = "dummy.trr"
    sample.init_diffusion_vacf("output/test/vacf.obj", len_correration=1e-14,
                               new_time_origin=2e-15, sample_step=1, len_frame=1e-15,
                               bin_num=bin_num, direction=direction)
    assert sample._is_diffusion_vacf
    return sample


def _run_vacf(sample, com):
    """Drive _diffusion_vacf once over a single frame of centres of mass."""
    corr_steps = sample._diff_vacf_inp["corr_steps"]
    com = np.asarray(com, float)
    pos_list = np.repeat(com[np.newaxis, :, :], corr_steps, axis=0)
    vel_list = np.ones((2*corr_steps-1, len(com), 3), float)
    data = sample._diffusion_vacf_data()
    sample._diffusion_vacf(data, corr_steps, pos_list, 0, vel_list, 0)
    return data["shape_00"]


def _pore_geometry(sample):
    """Focal axis, radius and mid-pore z of the single shape."""
    focal = sample._pore_props["shape_00"]["focal"]
    return focal, sample._pore_props["shape_00"]["diam"]/2, \
        focal[2] + sample._pore_props["box"]["res"]


def test_diffusion_vacf_radial_bins_by_distance(change_to_tests_dir):
    """Radial sampling must bin by the distance to the pore axis.

    The binning coordinate was pos_list[..., 0], the raw x coordinate of the
    box, so molecules were binned by x instead of by their distance from the
    cylinder axis.
    """
    sample = _vacf_sample(bin_num=4)
    assert sample._diff_vacf_inp["direction"] == "radial_cylindrical"
    focal, radius, z_mid = _pore_geometry(sample)

    # One molecule per radial bin; the second is offset along y, not x, so a
    # binning that reads the x coordinate cannot place it in bin 1
    com = [[focal[0] + (i + 0.5)/4*radius, focal[1], z_mid] for i in range(4)]
    com[1] = [focal[0], focal[1] + 1.5/4*radius, z_mid]
    assert list(_run_vacf(sample, com)["density"][:, 0]) == [1, 1, 1, 1]


def test_diffusion_vacf_wall_mask_uses_radial_distance(change_to_tests_dir):
    """The pore wall test must use the distance to the shape's focal axis.

    in_wall_mask compared the binning coordinate against diam/2. That
    coordinate was an absolute box coordinate, so a molecule well inside the
    pore was discarded whenever the pore itself sits far from the origin.
    """
    sample = _vacf_sample(bin_num=4)
    focal, radius, z_mid = _pore_geometry(sample)

    # Inside the pore, but at an x far larger than the pore radius
    assert _run_vacf(sample, [[focal[0] + radius*0.5, focal[1], z_mid]])["density"].sum() == 1
    # Beyond the wall
    assert _run_vacf(sample, [[focal[0] + radius*2.0, focal[1], z_mid]])["density"].sum() == 0


def test_diffusion_vacf_pore_cartesian_direction(change_to_tests_dir):
    """A pore system must also be sampleable along a box axis.

    _bin_box returns one bin list shared by every shape while _bin_pore
    returns one per shape, but the sampler indexed the bins by pore_id in
    both cases, so a cartesian direction raised
    TypeError: list indices must be integers or slices, not str.
    """
    sample = _vacf_sample(bin_num=4, direction=2)
    bins = sample._diff_vacf_inp["bins"]
    assert isinstance(bins, list)
    focal, radius, z_mid = _pore_geometry(sample)

    # Two molecules inside the pore, in the two z bins the pore spans, plus
    # one at the same z but out in the pore wall
    inside = focal[0] + radius*0.5
    density = _run_vacf(sample, [[inside, focal[1], 7.0],
                                 [inside, focal[1], 12.0],
                                 [focal[0] + radius*2.0, focal[1], 7.0]])["density"][:, 0]
    assert list(density) == [0, 1, 1, 0]


def _write_vacf_sample(link, bin_num=4):
    """Write a VACF sample file whose first bin is occupied and the rest empty."""
    sample = _vacf_sample(bin_num=bin_num)
    data = sample._diffusion_vacf_data()

    # Bin 0 saw two molecules with a constant correlation of 2, so the
    # normalised vacf of that bin is exactly 1 in every direction
    data["shape_00"]["vacf_data"][0, :, 0, :] = 2.0
    data["shape_00"]["density"][0, 0] = 2

    inp = dict(sample._diff_vacf_inp)
    inp["num_frame"] = sample._num_frame
    inp["mass"] = 1.0
    pa.utils.mkdirp("output/test")
    pa.utils.save({"type": "diff_vacf", "pore": sample._pore,
                   "inp": inp, "data": data}, link)
    return sample._diff_vacf_inp


def test_integrate_bin_diffusion_vacf_unoccupied_bins_are_zero(change_to_tests_dir):
    """An unoccupied bin must integrate to 0.0, not NaN.

    The normalisation divided by the per-bin density without guarding
    against an empty bin, so 0/0 put NaN in every direction of every bin no
    molecule ever entered, and the NaN then spread through any mean.
    """
    link = "output/test/vacf_zero.obj"
    inp = _write_vacf_sample(link, bin_num=4)
    integrated = pa.diffusion.integrate_bin_diffusion_vacf(link, pore_id="shape_00")

    # bins x residues x correlation steps x directions
    assert integrated.shape == (4, 1, inp["corr_steps"], 3)
    assert np.all(np.isfinite(integrated))

    # The three empty bins contribute nothing at all
    assert not np.any(integrated[1:])

    # The occupied bin integrates a constant 1 over the correlation time
    dt = inp["len_frame"]*inp["sample_step"]
    expect = np.arange(inp["corr_steps"])*dt
    for axis in range(3):
        assert np.allclose(integrated[0, 0, :, axis], expect)


def test_init_diffusion_vacf_rejects_bad_direction(change_to_tests_dir):
    """An unusable direction must be refused.

    The guard reads `direction not in [0,1,2] and (direction == "radial" and
    not self._pore)`, whose right operand is False for anything that is not
    "radial", so every other wrong value is accepted.
    """
    mol = pms.Molecule(inp="data/benzene.gro")
    for direction in (5, "axial"):
        sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
        sample._traj = "dummy.trr"
        sample.init_diffusion_vacf("output/test/vacf.obj", direction=direction)
        assert not sample._is_diffusion_vacf, direction
