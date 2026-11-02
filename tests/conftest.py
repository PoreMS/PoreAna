import os
import shutil

import pytest


@pytest.fixture(autouse=True, scope="session")
def change_to_tests_dir():
    """Run the suite from inside tests/ so the relative data paths resolve.

    Unlike the old setUpClass, this restores the original working directory on
    teardown, so coverage and pytest artefacts stay where they were invoked.
    """
    original = os.getcwd()
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(tests_dir)
    yield
    os.chdir(original)


@pytest.fixture(scope="session")
def sample_output(change_to_tests_dir):
    """Clear tests/output and run every trajectory sampling once per session.

    The unittest suite repeated this block for every test class; running it once
    is what makes the split into unit and integration tests worthwhile.
    """
    import porems as pms
    import poreana as pa

    folder = "output"
    pa.utils.mkdirp(folder)
    pa.utils.mkdirp(folder+"/temp")
    open(folder+"/temp.txt", "a").close()

    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)

    # Load molecules
    mol_B = pms.Molecule(inp="data/benzene.gro")
    mol_W = pms.Molecule(inp="data/spc216.gro")
    mol_H = pms.Molecule(inp="data/heptane.gro")

    len_step_mc = [1, 2, 5, 10, 20, 30, 40, 50, 100, 200, 250, 300, 350]

    # Single core
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_density("output/dens_cyl_s.obj")
    sample.init_gyration("output/gyr_cyl_s.obj")
    sample.init_angle("output/angle_cyl.obj", [0, 3])
    sample.sample(is_parallel=False)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_density("output/dens_cyl_s_const_a.obj", bin_const_A=True)
    sample.sample(is_parallel=False)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_diffusion_bin("output/diff_cyl_s.obj")
    sample.sample(is_parallel=False)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_density("output/dens_cyl_no_remove.obj", remove_pore_from_res=False)
    sample.sample(is_parallel=False)

    sample = pa.Sample("data/pore_system_slit_new.yml", "data/traj_slit.xtc", mol_W)
    sample.init_density("output/dens_slit.obj")
    sample.sample(is_parallel=False, is_pbc=False)

    sample = pa.Sample([6.00035, 6.00035, 19.09191], "data/traj_box.xtc", mol_H)
    sample.init_density("output/dens_box.obj")
    sample.init_gyration("output/gyr_box.obj")
    sample.init_angle("output/angle_box.obj", [0, 4])
    sample.sample(shift=[0, 0, 3.3], is_parallel=False)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_diffusion_mc("output/diff_mc_cyl_s.obj", len_step=len_step_mc)
    sample.sample(is_parallel=False)

    sample = pa.Sample([6.00035, 6.00035, 19.09191], "data/traj_box.xtc", mol_H)
    sample.init_diffusion_mc("output/diff_mc_box.obj", len_step=[1, 2, 5, 10, 20, 30, 40, 50])
    sample.sample(shift=[0, 0, 3.3], is_parallel=False, is_pbc=True)

    # Parallel
    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_density("output/dens_cyl_p.obj")
    sample.init_gyration("output/gyr_cyl_p.obj")
    sample.init_diffusion_bin("output/diff_cyl_p.obj")
    sample.sample(is_parallel=True, is_pbc=False)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_diffusion_mc("output/diff_mc_cyl_p.obj", len_step=len_step_mc)
    sample.sample(is_parallel=True, is_pbc=True, num_cores=6)

    # VACF needs velocities, so a trr. traj_oxy.trr holds 40 O2 molecules in a
    # two phase box, 10001 frames written every 20 fs over 200 ps. The box is
    # read from the trajectory, the lower quarter of it holds no oxygen at all.
    mol_O = pms.Molecule("oxygen", "OXY")
    mol_O.add("O", [0, 0, 0], name="O1")
    mol_O.add("O", [0, 0, 0.1208], name="O2")

    # len_frame is the md timestep and sample_step the md steps per written
    # frame, so a frame is 20 fs apart and the correlation window is 4 ps
    vacf_inp = {"len_correration": 4e-12, "new_time_origin": 4e-13,
                "sample_step": 20, "len_frame": 1e-15, "bin_num": 20, "direction": 2}

    sample = pa.Sample([], "data/traj_oxy.trr", mol_O)
    sample.init_diffusion_vacf("output/vacf_oxy.obj", **vacf_inp)
    sample.sample(is_parallel=False, is_pbc=True)

    sample = pa.Sample([], "data/traj_oxy.trr", mol_O)
    sample.init_diffusion_vacf("output/vacf_oxy_p.obj", **vacf_inp)
    sample.sample(is_parallel=True, is_pbc=True)

    sample = pa.Sample([], "data/traj_oxy.trr", mol_O)
    sample.init_diffusion_vacf("output/vacf_oxy_res.obj", sample_each_residue=True, **vacf_inp)
    sample.sample(is_parallel=False, is_pbc=True)

    # Small MC output for the file_to_text and table tests. Neither asserts a
    # physical value, so the step counts only have to produce a loadable file.
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    model._len_step = [10, 20, 40]
    pa.MC().run(model, "output/mc_fixture.obj", nmc_eq=20, nmc=20, is_print=False, is_parallel=False)

    return True
