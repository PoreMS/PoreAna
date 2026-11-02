"""Tests that run against the trajectories sampled by the session fixture."""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

import porems as pms
import poreana as pa


##########
# Sample #
##########
def test_sample(sample_output):
    mol = pms.Molecule("benzene", "BEN", inp="data/benzene.gro")
    mol2 = pms.Molecule("benzene", "BEN", inp="data/benzene.gro")
    mol2.add("H", [0, 0, 0])
    mol3 = pms.Molecule("water", "SOL", inp="data/spc216.gro")

    # Inconsistent atom and mass counts, and a wrong shift dimension
    pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol2)
    pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol, atoms=["C1"], masses=[1, 1])
    pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol).sample(shift=[1])

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol, atoms=["C1"])
    sample.init_diffusion_bin("output/diff_np_s.obj", len_obs=3e-12)

    sample = pa.Sample([0, 0, 1], "data/traj_cylinder.xtc", mol, atoms=["C1"])
    sample.init_diffusion_bin("output/diff_box_test.obj", len_obs=3e-12)

    # Bin and MC diffusion exclude each other, in both orders
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_diffusion_bin("output/test.obj")
    sample.init_diffusion_mc("output/test.obj", len_step=[1, 2, 5, 10])
    assert not sample._is_diffusion_mc

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_diffusion_mc("output/test.obj", len_step=[1, 2, 5, 10])
    sample.init_diffusion_bin("output/test.obj")
    assert not sample._is_diffusion_bin

    # Wrong direction
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_diffusion_mc("output/test.obj", len_step=[1, 2, 5, 10], direction=4)
    assert not sample._is_diffusion_mc

    # Angle cannot be parallelized
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_angle("output/angle_cyl.obj", [0, 3])
    sample.sample(is_parallel=True, is_pbc=False)

    # Shape without a predefined normal
    sample = pa.Sample("data/pore_system_slit_new.yml", "data/traj_slit.xtc", mol3)
    sample.init_angle("output/angle_cyl.obj", [0, 3])


def test_sample_partial_atom_selection(sample_output):
    """A subset of a molecule's atoms must be sampleable.

    _sample_helper reshapes the whole frame to (num_res, len(self._atoms), 3),
    but a frame holds every atom of every residue, so any partial selection
    gives a reshape of the wrong size.
    """
    mol = pms.Molecule("benzene", "BEN", inp="data/benzene.gro")
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc",
                       mol, atoms=["C1"], masses=[1])
    sample.init_density("output/dens_partial.obj")
    sample.sample(is_parallel=False)
    assert os.path.isfile("output/dens_partial.obj")


def test_sample_numpy_positions_only(sample_output):
    """Sampling positions alone must not touch the velocities.

    _sample_helper read frame.velocities for every numpy run regardless of
    the velocities flag, so a positions-only run against an xtc raised
    ChemfilesError: velocity data is not defined in this frame.
    """
    mol = pms.Molecule("benzene", "BEN", inp="data/benzene.gro")
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_numpy_file("output/numpy_pos.npz", positions=True, velocities=False)
    sample.sample(is_parallel=False, is_pbc=False)

    data = np.load("output/numpy_pos.npz")
    assert list(data.keys()) == ["positions"]
    assert data["positions"].shape == (sample._num_frame, sample.num_res, 3)


###########
# Density #
###########
def test_density(sample_output):
    dens_s = pa.density.bins("output/dens_cyl_s.obj", target_dens=16)
    dens_p = pa.density.bins("output/dens_cyl_p.obj", target_dens=16)
    dens_no_remove = pa.density.bins("output/dens_cyl_no_remove.obj")
    dens_slit = pa.density.bins("output/dens_slit.obj", target_dens=997)
    dens_box = pa.density.bins("output/dens_box.obj")
    pa.density.bins("output/dens_cyl_s_const_a.obj")

    plt.figure()
    pa.density.bins_plot(dens_s, target_dens=0.146, pore_id="shape_00")
    pa.density.bins_plot(dens_s, target_dens=0.146, pore_id="shape_00", is_mean=True)
    plt.savefig("output/density.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.density.bins_plot(dens_no_remove, pore_id="shape_00")
    plt.savefig("output/density_no_remove.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.density.bins_plot(dens_slit, is_mean=True, pore_id="shape_00")
    plt.savefig("output/density_slit.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.density.bins_plot(dens_s, intent="in", pore_id="shape_00")
    pa.density.bins_plot(dens_p, intent="in", pore_id="shape_00")
    plt.savefig("output/density_in.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.density.bins_plot(dens_s, intent="ex", pore_id="shape_00")
    pa.density.bins_plot(dens_p, intent="ex", pore_id="shape_00")
    plt.savefig("output/density_ex.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.density.bins_plot(dens_box, intent="ex", pore_id="shape_00")
    plt.savefig("output/density_box.pdf", format="pdf", dpi=100)
    assert pa.density.bins_plot(dens_s, intent="DOTA") is None
    plt.close("all")

    assert round(dens_s["dens"]["shape_00"]["in"], 3) == 13.01
    assert round(dens_s["dens"]["ex"], 3) == 15.707
    assert round(dens_p["dens"]["shape_00"]["in"], 3) == 13.01
    assert round(dens_p["dens"]["ex"], 3) == 15.707

    # Serial and parallel sampling must agree
    assert dens_s["num_dens"]["shape_00"]["in"] == dens_p["num_dens"]["shape_00"]["in"]


def test_density_mean(sample_output):
    dens_s = pa.density.bins("output/dens_cyl_s.obj", is_print=False)
    mean = pa.density.mean(dens_s)
    assert set(mean) == {"num_dens", "dens"}
    assert round(mean["num_dens"]["shape_00"], 6) == 0.127489
    assert round(mean["dens"]["shape_00"], 3) == 16.537


##############
# Adsorption #
##############
def test_adsorption(sample_output):
    ads_s = pa.adsorption.calculate("output/dens_cyl_s.obj")
    ads_p = pa.adsorption.calculate("output/dens_cyl_p.obj")

    assert round(ads_s["conc"]["shape_00"]["mumol_m2"], 4) == 0.1514
    assert round(ads_s["conc"]["mmol_l"], 2) == 212.97
    assert round(ads_s["num"]["shape_00"]["in"], 4) == 10.5912
    assert round(ads_s["num"]["ex"], 4) == 28.6817
    assert round(ads_p["num"]["shape_00"]["in"], 4) == 10.5912


def test_adsorption_does_not_mutate_the_loaded_sample(sample_output, monkeypatch):
    """Repeated calculations over one loaded sample must agree.

    calculate subtracted the two reservoirs from pore["box"]["dimensions"]
    in place. It reloads the file on every call, so the damage was invisible
    from the outside, but any caller or loader that hands it the same
    dictionary twice saw the box shrink on each call.
    """
    sample = pa.utils.load("output/dens_cyl_s.obj")
    box_before = list(sample["pore"]["box"]["dimensions"])
    monkeypatch.setattr(pa.adsorption.utils, "load", lambda *args, **kwargs: sample)

    first = pa.adsorption.calculate("output/dens_cyl_s.obj")
    second = pa.adsorption.calculate("output/dens_cyl_s.obj")

    assert list(sample["pore"]["box"]["dimensions"]) == box_before
    assert first["conc"]["shape_00"]["mumol_m2"] == second["conc"]["shape_00"]["mumol_m2"]
    assert round(first["conc"]["shape_00"]["mumol_m2"], 4) == 0.1514


############
# Gyration #
############
def test_gyration(sample_output):
    plt.figure()
    mean_s = pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", is_mean=True)
    plt.savefig("output/gyration_s.pdf", format="pdf", dpi=100)
    plt.figure()
    mean_p = pa.gyration.bins_plot("output/gyr_cyl_p.obj", "output/dens_cyl_p.obj")
    plt.savefig("output/gyration_p.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="in")
    pa.gyration.bins_plot("output/gyr_cyl_p.obj", "output/dens_cyl_p.obj", intent="in")
    plt.savefig("output/gyration_in.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="ex")
    pa.gyration.bins_plot("output/gyr_box.obj", "output/dens_box.obj", intent="ex")
    plt.savefig("output/gyration_ex.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="in", is_norm=True)
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="ex", is_norm=True)
    plt.savefig("output/gyration_norm.pdf", format="pdf", dpi=100)
    assert pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="DOTA") is None
    plt.close("all")

    assert round(mean_s["in"], 2) == 0.13
    assert round(mean_s["ex"], 2) == 0.15
    assert round(mean_p["in"], 2) == 0.13
    assert round(mean_p["ex"], 2) == 0.15


#########
# Angle #
#########
def test_angle(sample_output):
    plt.figure()
    mean_a = pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", is_mean=True)
    plt.savefig("output/angle.pdf", format="pdf", dpi=100)
    plt.figure()
    pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", intent="in", is_norm=True)
    pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", intent="ex", is_norm=True)
    plt.savefig("output/angle_norm.pdf", format="pdf", dpi=100)
    assert pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", intent="DOTA") is None
    plt.close("all")

    assert "in" in mean_a
    assert "ex" in mean_a


def test_angle_gyration_box_system(sample_output):
    """angle and gyration must plot a box system, which has no interior.

    width["in"] is an empty list for a box system, and both functions index
    it; angle also reads data[pore_id]["in"] before checking for a pore.
    """
    plt.figure()
    mean_g = pa.gyration.bins_plot("output/gyr_box.obj", "output/dens_box.obj")
    mean_a = pa.angle.bins_plot("output/angle_box.obj", "output/dens_box.obj")
    plt.close("all")

    assert set(mean_g) == {"ex"}
    assert set(mean_a) == {"ex"}

    plt.figure()
    pa.gyration.bins_plot("output/gyr_box.obj", "output/dens_box.obj", is_norm=True)
    pa.angle.bins_plot("output/angle_box.obj", "output/dens_box.obj", is_norm=True)
    plt.close("all")


#################
# Bin diffusion #
#################
def test_diffusion_bin(sample_output):
    plt.figure()
    bins_s = pa.diffusion.bins("output/diff_cyl_s.obj")
    pa.diffusion.bins_plot(bins_s)
    pa.diffusion.bins_plot(pa.diffusion.bins("output/diff_cyl_s.obj", is_norm=True))
    plt.savefig("output/diffusion_bins.pdf", format="pdf", dpi=100)
    plt.close("all")

    assert set(bins_s) >= {"diff", "width", "is_norm"}
    assert "shape_00" in bins_s["diff"]
    assert len(bins_s["diff"]["shape_00"]) > 0

    plt.figure()
    mean_s = pa.diffusion.mean(pa.diffusion.bins("output/diff_cyl_s.obj"),
                               pa.density.bins("output/dens_cyl_s.obj"))
    mean_p = pa.diffusion.mean(pa.diffusion.bins("output/diff_cyl_p.obj"),
                               pa.density.bins("output/dens_cyl_p.obj"))
    plt.savefig("output/diff_mean_check.pdf", format="pdf", dpi=100)
    plt.close("all")
    assert round(mean_s["shape_00"], 2) == 1.13
    assert round(mean_p["shape_00"], 2) == 1.13


def test_diffusion_cui_axial(sample_output):
    plt.figure()
    pa.diffusion.cui("output/diff_cyl_s.obj", intent="axial", is_fit=True)
    plt.savefig("output/diffusion_cui.pdf", format="pdf", dpi=100)
    plt.close("all")


def test_diffusion_cui_radial(sample_output):
    """The radial fit must read the nested pore data model.

    It takes the starting radius from pore["diam"], which moved to
    pore[shape_id]["diam"], so the tutorial example in
    docs/source/diffusion_bin.md raises KeyError: 'diam'.
    """
    plt.figure()
    pa.diffusion.cui("output/diff_cyl_s.obj", is_fit=True)
    plt.close("all")


################
# MC diffusion #
################
def test_diffusion_mc_model(sample_output):
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    assert np.array_equal(np.round(model._diff_bin, 3), np.array([-3.702] * model._bin_num))
    assert np.array_equal(model._df_bin, np.array([0] * model._bin_num))

    model = pa.StepModel("output/diff_mc_cyl_s.obj", 6, 10)
    assert np.array_equal(np.round(model._diff_bin, 3), np.array([-3.702] * model._bin_num))
    assert np.array_equal(model._df_bin, np.array([0] * model._bin_num))


def test_diffusion_mc_likelihood(sample_output):
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    for len_step, expected in ((1, -149411.12), (2, -169616.56), (10, -234093.97)):
        pa.MC._len_step = len_step
        assert round(pa.MC()._log_likelihood_z(model), 2) == expected


# nmc_eq=1000 is the shortest equilibration that keeps the 0.3 tolerance
# meaningful. Over five repeats the fitted coefficient strays by at most 0.1;
# halving the steps again triples that spread.
MC_EQ = 1000
MC_PROD = 2000


def test_diffusion_mc_run_serial(sample_output):
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    model._len_step = [10, 20, 40]
    pa.MC().run(model, "output/diff_test_mc.yml", nmc_eq=MC_EQ, nmc=MC_PROD,
                is_print=False, is_parallel=False)

    plt.figure()
    diff = pa.diffusion.mc_fit("output/diff_test_mc.yml")
    plt.savefig("output/mc_fit.pdf", format="pdf", dpi=100)
    plt.figure()
    diff_pore = pa.diffusion.mc_fit("output/diff_test_mc.yml", section="pore")
    plt.savefig("output/mc_fit_pore.pdf", format="pdf", dpi=100)
    plt.close("all")

    assert abs(diff[0] - 1.6) < 0.3
    assert abs(diff_pore[0] - 1.2) < 0.3


def test_diffusion_mc_run_parallel(sample_output):
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    model._len_step = [10, 20, 40]
    pa.MC().run(model, "output/diff_test_mc.obj", nmc_eq=MC_EQ, nmc=MC_PROD,
                is_print=False, is_parallel=True)

    plt.figure()
    diff = pa.diffusion.mc_fit("output/diff_test_mc.obj")
    plt.figure()
    diff_pore = pa.diffusion.mc_fit("output/diff_test_mc.obj", section="pore")
    plt.close("all")

    assert abs(diff[0] - 1.6) < 0.3
    assert abs(diff_pore[0] - 1.2) < 0.3


def test_diffusion_mc_print_output(sample_output):
    """Cover the console output of both models and of a single lag time."""
    pa.StepModel("output/diff_mc_cyl_s.obj", 6, 10, is_print=True)
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10, is_print=True)
    model._len_step = [10]
    pa.MC().run(model, "output/diff_test_mc_spline.obj", nmc_eq=50, nmc=100,
                is_print=True, is_parallel=False)


def test_diffusion_mc_box(sample_output):
    model = pa.CosineModel("output/diff_mc_box.obj", 6, 10)
    model._len_step = [10, 20, 30, 40, 50]
    pa.MC().run(model, "output/diff_test_mc_box.obj", nmc_eq=100, nmc=300, is_print=False, is_parallel=False)
    plt.figure()
    diff = pa.diffusion.mc_fit("output/diff_test_mc_box.obj")
    plt.savefig("output/diffusion_fit_box.pdf", format="pdf", dpi=100)
    plt.close("all")
    assert abs(diff[0] - 1.3) < 0.3


def test_diffusion_mc_returns_coefficients(sample_output):
    """list_diff_coeff / list_df_coeff must hold the model coefficients.

    The output dictionary assigns list_diff_profile to both keys, so the
    per-step coefficient dictionaries that were built are thrown away.
    """
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    model._len_step = [10]
    pa.MC().run(model, "output/mc_coeff.obj", nmc_eq=100, nmc=100, is_print=False, is_parallel=False)

    out = pa.utils.load("output/mc_coeff.obj")["output"]
    # 6 cosine terms for the diffusion profile, 10 for the free energy
    assert len(out["list_diff_coeff"][10]) == 6
    assert len(out["list_df_coeff"][10]) == 10
    assert len(out["diff_profile"][10]) == model._bin_num


def test_diffusion_mc_parallel_fluctuations(sample_output):
    """Parallel MC must collect fluc_diff_bin / fluc_df_bin from every worker.

    The merge loop copies ten output keys from the workers but not these two,
    so the result only carries the lag times that landed on worker 0.
    """
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    model._len_step = [10, 20, 40]
    pa.MC().run(model, "output/mc_parallel.obj", nmc_eq=200, nmc=400,
                is_print=False, is_parallel=True, np=3)

    out = pa.utils.load("output/mc_parallel.obj")["output"]
    assert set(out["fluc_diff_bin"]) == {10, 20, 40}
    assert set(out["fluc_df_bin"]) == {10, 20, 40}


#####################
# Transition matrix #
#####################
def test_parallel_sample(sample_output):
    """Splitting the trajectory over workers must not change the result."""
    data_s = pa.utils.load("output/diff_mc_cyl_s.obj")["data"]
    data_p = pa.utils.load("output/diff_mc_cyl_p.obj")["data"]

    assert set(data_s) == set(data_p)
    for step in data_s:
        assert np.array_equal(data_s[step], data_p[step]), step

    # Anchor both against a committed serial reference, so a change that moves
    # serial and parallel together is caught as well
    reference = pa.utils.load("data/check_output_sample.obj")
    pore = pa.utils.load("data/pore_system_cylinder_new.yml", file_type="yml")
    assert reference["inp"]["bins"][-1] == pore["system"]["dimensions"][2]
    for step in reference["data"]:
        assert np.array_equal(reference["data"][step], data_s[step]), step
        assert np.array_equal(reference["data"][step], data_p[step]), step


###############
# Free energy #
###############
def test_freeenergy_mc(sample_output):
    plt.figure()
    pa.freeenergy.mc_profile("output/diff_test_mc.obj", is_plot=True)
    plt.savefig("output/energy_profile.pdf", format="pdf", dpi=100)
    plt.close("all")


##########
# Tables #
##########
def test_tables(sample_output):
    assert pa.tables.mc_model("output/mc_fixture.obj", print_con=True) is not None
    assert pa.tables.mc_inputs("output/mc_fixture.obj", print_con=True) is not None


def test_tables_mc_results(sample_output):
    """mc_results must read the nested pore data model.

    It still indexes data["pore"]["box"][2] and data["pore"]["res"], but box
    became a dictionary and res moved inside it.
    """
    results = pa.tables.mc_results("output/mc_fixture.obj", print_con=True)
    assert "pore" in results.index
    assert pa.tables.mc_results("output/mc_fixture.obj", sections={"DOTA": []}) is None


################
# Text export  #
################
def test_file_to_text(sample_output):
    """All four export branches must read the nested pore data model.

    Each one still reads pore["res"], pore["diam"], pore["box"] as a list and
    pore["type"], so every branch raises KeyError: 'res'.
    """
    pa.utils.file_to_text("output/dens_cyl_s.obj", "output/dens_cyl_s.txt")
    with open("output/dens_cyl_s.txt") as f:
        content = f.read()
    assert "[System]" in content
    assert "[Density]" in content
    assert "[Adsorption]" in content

    pa.utils.file_to_text("output/dens_box.obj", "output/dens_box.txt")
    with open("output/dens_box.txt") as f:
        content_box = f.read()
    assert "[Density]" in content_box
    assert "[Adsorption]" not in content_box

    pa.utils.file_to_text("output/diff_cyl_s.obj", "output/diff_cyl_s.txt",
                          link_dens="output/dens_cyl_s.obj")
    with open("output/diff_cyl_s.txt") as f:
        assert "[Diffusion]" in f.read()

    # A missing density link is reported, not raised
    pa.utils.file_to_text("output/diff_cyl_s.obj", "output/diff_missing_dens.txt")
    assert not os.path.isfile("output/diff_missing_dens.txt")

    pa.utils.file_to_text("output/gyr_cyl_s.obj", "output/gyr_cyl_s.txt",
                          link_dens="output/dens_cyl_s.obj")
    with open("output/gyr_cyl_s.txt") as f:
        assert "[Gyration]" in f.read()

    pa.utils.file_to_text("output/mc_fixture.obj", "output/mc_fixture.txt")
    with open("output/mc_fixture.txt") as f:
        assert "[Profiles]" in f.read()
