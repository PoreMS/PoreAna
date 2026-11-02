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


@pytest.mark.xfail(strict=True, reason="acos is called on a cosine that can exceed 1")
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
@pytest.mark.xfail(strict=True, reason="the outermost radial bin is dropped")
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
@pytest.mark.xfail(strict=True, reason="init_diffusion_mc sorts the caller's list")
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


@pytest.mark.xfail(strict=True, reason="type(d0) == float is False for np.float64")
def test_model_accepts_numpy_float():
    """A numpy float guess must initialize the diffusion profile.

    The guard is `type(d0) == float`, which no numpy scalar satisfies, so
    self._d0 is never assigned at all.
    """
    link = _write_minimal_mc_sample("output/test/mc_minimal.obj")
    model = pa.CosineModel(link, 2, 2, d0=np.float64(1.0))
    assert np.array_equal(np.round(model._diff_bin, 4), np.array([-6.9078] * model._bin_num))


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


@pytest.mark.xfail(strict=True, reason="the direction guard can never be true")
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
