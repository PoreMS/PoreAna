import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import porems as pms
import poreana as pa


# ---------------------------------------------------------------------------
# Sample
# ---------------------------------------------------------------------------

def test_sample(sample_output):
    mol = pms.Molecule("benzene", "BEN", inp="data/benzene.gro")
    mol2 = pms.Molecule("benzene", "BEN", inp="data/benzene.gro")
    mol2.add("H", [0, 0, 0])
    mol3 = pms.Molecule("water", "SOL", inp="data/spc216.gro")

    pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol2)
    pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol, atoms=["C1"], masses=[1, 1])
    pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol).sample(shift=[1])

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol, atoms=["C1"])
    sample.init_diffusion_bin("output/diff_np_s.obj", len_obs=3e-12)

    sample = pa.Sample([0, 0, 1], "data/traj_cylinder.xtc", mol, atoms=["C1"])
    sample.init_diffusion_bin("output/diff_box_test.obj", len_obs=3e-12)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_diffusion_bin("output/test.obj")
    sample.init_diffusion_mc("output/test.obj", len_step=[1, 2, 5, 10])

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_diffusion_mc("output/test.obj", len_step=[1, 2, 5, 10])
    sample.init_diffusion_bin("output/test.obj")

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_diffusion_mc("output/test.obj", len_step=[1, 2, 5, 10], direction=4)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol)
    sample.init_angle("output/angle_cyl.obj", [0, 3])
    sample.sample(is_parallel=True, is_pbc=False)

    sample = pa.Sample("data/pore_system_slit_new.yml", "data/traj_slit.xtc", mol3)
    sample.init_angle("output/angle_cyl.obj", [0, 3])


# ---------------------------------------------------------------------------
# Density
# ---------------------------------------------------------------------------

def test_density(sample_output):
    dens_s = pa.density.bins("output/dens_cyl_s.obj", target_dens=16)
    dens_p = pa.density.bins("output/dens_cyl_p.obj", target_dens=16)
    mean_w = pa.density.mean(dens_s)
    pa.density.bins("output/dens_cyl_no_remove.obj")
    pa.density.bins("output/dens_slit.obj", target_dens=997)
    pa.density.bins("output/dens_box.obj")

    plt.figure()
    pa.density.bins_plot(dens_s, target_dens=0.146, pore_id="shape_00")
    pa.density.bins_plot(dens_s, target_dens=0.146, pore_id="shape_00", is_mean=True)
    plt.savefig("output/density.pdf", format="pdf", dpi=100)
    plt.close("all")

    plt.figure()
    pa.density.bins_plot(dens_s, intent="in", pore_id="shape_00")
    pa.density.bins_plot(dens_s, intent="ex")
    pa.density.bins_plot(dens_s, intent="DOTA")
    plt.close("all")

    # adsorption
    ads = pa.adsorption.calculate("output/dens_cyl_s.obj")
    assert "conc" in ads
    assert "num" in ads
    assert "shape_00" in ads["conc"]
    assert "mumol_m2" in ads["conc"]["shape_00"]
    assert "mmol_l" in ads["conc"]
    assert ads["num"]["ex"] > 0

    # density.mean weighted integration
    assert "num_dens" in mean_w
    assert "dens" in mean_w
    assert "shape_00" in mean_w["num_dens"]
    assert mean_w["num_dens"]["shape_00"] > 0

    assert round(dens_s["dens"]["shape_00"]["in"], 3) == 13.01
    assert round(dens_s["dens"]["ex"], 3) == 15.707
    assert round(dens_p["dens"]["shape_00"]["in"], 3) == 13.01
    assert round(dens_p["dens"]["ex"], 3) == 15.707


# ---------------------------------------------------------------------------
# Gyration
# ---------------------------------------------------------------------------

def test_gyration(sample_output):
    plt.figure()
    mean_s = pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", is_mean=True)
    plt.savefig("output/gyration_s.pdf", format="pdf", dpi=100)
    plt.close("all")

    plt.figure()
    mean_p = pa.gyration.bins_plot("output/gyr_cyl_p.obj", "output/dens_cyl_p.obj")
    plt.savefig("output/gyration_p.pdf", format="pdf", dpi=100)
    plt.close("all")

    plt.figure()
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="in")
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="ex")
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="DOTA")
    plt.close("all")

    # normalized x-axis
    plt.figure()
    pa.gyration.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", is_norm=True)
    plt.close("all")

    assert round(mean_s["in"], 2) == 0.13
    assert round(mean_s["ex"], 2) == 0.15
    assert round(mean_p["in"], 2) == 0.13
    assert round(mean_p["ex"], 2) == 0.15


# ---------------------------------------------------------------------------
# Angle
# ---------------------------------------------------------------------------

def test_angle(sample_output):
    plt.figure()
    mean_a = pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", is_mean=True)
    plt.savefig("output/angle.pdf", format="pdf", dpi=100)
    plt.close("all")

    plt.figure()
    pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", intent="in")
    pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", intent="ex")
    pa.angle.bins_plot("output/gyr_cyl_s.obj", "output/dens_cyl_s.obj", intent="DOTA")
    plt.close("all")

    # normalized x-axis
    plt.figure()
    pa.angle.bins_plot("output/angle_cyl.obj", "output/dens_cyl_s.obj", is_norm=True)
    plt.close("all")

    assert "in" in mean_a
    assert "ex" in mean_a


# ---------------------------------------------------------------------------
# Bin diffusion
# ---------------------------------------------------------------------------

def test_diffusion_bin(sample_output):
    plt.figure()
    bins_s = pa.diffusion.bins("output/diff_cyl_s.obj")
    pa.diffusion.bins_plot(bins_s)
    pa.diffusion.bins_plot(pa.diffusion.bins("output/diff_cyl_s.obj", is_norm=True))
    plt.savefig("output/diffusion_bins.pdf", format="pdf", dpi=100)
    plt.close("all")

    # output structure
    assert "diff" in bins_s
    assert "width" in bins_s
    assert "is_norm" in bins_s
    assert isinstance(bins_s["diff"], dict)
    assert "shape_00" in bins_s["diff"]
    assert len(bins_s["diff"]["shape_00"]) > 0

    mean_s = pa.diffusion.mean(
        pa.diffusion.bins("output/diff_cyl_s.obj"),
        pa.density.bins("output/dens_cyl_s.obj"),
    )
    mean_p = pa.diffusion.mean(
        pa.diffusion.bins("output/diff_cyl_p.obj"),
        pa.density.bins("output/dens_cyl_p.obj"),
    )
    assert round(mean_s["shape_00"], 2) == 1.13
    assert round(mean_p["shape_00"], 2) == 1.13


# ---------------------------------------------------------------------------
# MC diffusion — model
# ---------------------------------------------------------------------------

def test_diffusion_mc_model(sample_output):
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    assert np.array_equal(np.round(model._diff_bin, 3), np.array([-3.702] * model._bin_num))
    assert np.array_equal(model._df_bin, np.array([0] * model._bin_num))

    model = pa.StepModel("output/diff_mc_cyl_s.obj", 6, 10)
    assert np.array_equal(np.round(model._diff_bin, 3), np.array([-3.702] * model._bin_num))
    assert np.array_equal(model._df_bin, np.array([0] * model._bin_num))


# ---------------------------------------------------------------------------
# MC diffusion — MC run
# ---------------------------------------------------------------------------

def test_diffusion_mc_mc(sample_output):
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    pa.MC._len_step = 1
    assert round(pa.MC()._log_likelihood_z(model), 2) == -149411.12
    pa.MC._len_step = 2
    assert round(pa.MC()._log_likelihood_z(model), 2) == -169616.56
    pa.MC._len_step = 10
    assert round(pa.MC()._log_likelihood_z(model), 2) == -234093.97

    model._len_step = [10, 20, 40]
    pa.MC().run(model, "output/diff_test_mc.yml", nmc_eq=1000, nmc=2000, is_print=False, is_parallel=False)

    # verify list_diff_coeff holds model coefficients (one per cosine term),
    # not per-bin profile values — this was the bug that was fixed
    mc_out = pa.utils.load("output/diff_test_mc.yml")
    out = mc_out["output"]
    assert "list_diff_coeff" in out
    assert "diff_profile" in out
    step = model._len_step[0]
    # CosineModel with 6 terms: coeff list has 6 elements, profile has bin_num elements
    coeff_list = out["list_diff_coeff"][step]
    profile_list = out["diff_profile"][step]
    assert len(coeff_list) == 6
    assert len(profile_list) != 6

    diff = pa.diffusion.mc_fit("output/diff_test_mc.yml")
    diff_pore = pa.diffusion.mc_fit("output/diff_test_mc.yml", section="pore")
    assert abs(diff[0] - 1.6) < 0.3
    assert abs(diff_pore[0] - 1.2) < 0.3
    plt.close("all")

    pa.MC().run(model, "output/diff_test_mc.obj", nmc_eq=8000, nmc=2000, is_print=False, is_parallel=True)
    diff = pa.diffusion.mc_fit("output/diff_test_mc.obj")
    diff_pore = pa.diffusion.mc_fit("output/diff_test_mc.obj", section="pore")
    assert abs(diff[0] - 1.6) < 0.3
    assert abs(diff_pore[0] - 1.2) < 0.3
    plt.close("all")


# ---------------------------------------------------------------------------
# MC diffusion — box
# ---------------------------------------------------------------------------

def test_diffusion_mc_box(sample_output):
    model = pa.CosineModel("output/diff_mc_box.obj", 6, 10)
    model._len_step = [10, 20, 30, 40, 50]
    pa.MC().run(model, "output/diff_test_mc_box.obj", nmc_eq=100, nmc=300, is_print=False, is_parallel=False)
    diff = pa.diffusion.mc_fit("output/diff_test_mc_box.obj")
    assert abs(diff[0] - 1.3) < 0.3
    plt.close("all")


# ---------------------------------------------------------------------------
# Parallel sample consistency
# ---------------------------------------------------------------------------

def test_parallel_sample(sample_output):
    data_s = pa.utils.load("output/diff_mc_cyl_s.obj")["data"]
    data_p = pa.utils.load("output/diff_mc_cyl_p.obj")["data"]
    data_check = pa.utils.load("data/check_output_sample.obj")["data"]
    for i in [1, 2, 5, 10, 20, 30, 40]:
        assert np.array_equal(data_s[i], data_p[i])
        assert np.array_equal(data_check[i], data_p[i])


# ---------------------------------------------------------------------------
# file_to_text export
# ---------------------------------------------------------------------------

def test_file_to_text(sample_output):
    import os

    # dens_bin — pore
    pa.utils.file_to_text("output/dens_cyl_s.obj", "output/dens_cyl_s.txt")
    assert os.path.isfile("output/dens_cyl_s.txt")
    with open("output/dens_cyl_s.txt") as f:
        content = f.read()
    assert "[System]" in content
    assert "[Density]" in content
    assert "[Adsorption]" in content

    # dens_bin — box (no pore, no adsorption section)
    pa.utils.file_to_text("output/dens_box.obj", "output/dens_box.txt")
    assert os.path.isfile("output/dens_box.txt")
    with open("output/dens_box.txt") as f:
        content_box = f.read()
    assert "[Density]" in content_box
    assert "[Adsorption]" not in content_box

    # diff_bin — requires density link
    pa.utils.file_to_text("output/diff_cyl_s.obj", "output/diff_cyl_s.txt",
                          link_dens="output/dens_cyl_s.obj")
    assert os.path.isfile("output/diff_cyl_s.txt")
    with open("output/diff_cyl_s.txt") as f:
        assert "[Diffusion]" in f.read()

    # diff_bin — missing dens link exits gracefully (no file written, no crash)
    pa.utils.file_to_text("output/diff_cyl_s.obj", "output/diff_missing_dens.txt")

    # gyr_bin — requires density link
    pa.utils.file_to_text("output/gyr_cyl_s.obj", "output/gyr_cyl_s.txt",
                          link_dens="output/dens_cyl_s.obj")
    assert os.path.isfile("output/gyr_cyl_s.txt")
    with open("output/gyr_cyl_s.txt") as f:
        assert "[Gyration]" in f.read()

    # mc — full MC output (uses fixture-created mc_fixture.yml)
    pa.utils.file_to_text("output/mc_fixture.yml", "output/mc_fixture.txt")
    assert os.path.isfile("output/mc_fixture.txt")
    with open("output/mc_fixture.txt") as f:
        content_mc = f.read()
    assert "[Profiles]" in content_mc


# ---------------------------------------------------------------------------
# Free energy
# ---------------------------------------------------------------------------

def test_freeenergy_mc(sample_output):
    plt.figure()
    pa.freeenergy.mc_profile("data/check_output.obj")
    plt.savefig("output/energy_profile.pdf", format="pdf", dpi=100)
    plt.close("all")
