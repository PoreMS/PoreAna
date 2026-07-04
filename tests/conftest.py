import os
import shutil

import pytest


@pytest.fixture(autouse=True, scope="session")
def change_to_tests_dir():
    """Run all tests from within the tests/ directory so relative paths work."""
    original = os.getcwd()
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(tests_dir)
    yield
    os.chdir(original)


@pytest.fixture(scope="session")
def sample_output(change_to_tests_dir):
    """Session fixture: clean output dir and run all trajectory sampling once."""
    import porems as pms
    import poreana as pa

    folder = "output"
    pa.utils.mkdirp(folder)
    pa.utils.mkdirp(folder + "/temp")
    open(folder + "/temp.txt", "a").close()
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)

    mol_B = pms.Molecule(inp="data/benzene.gro")
    mol_W = pms.Molecule(inp="data/spc216.gro")
    mol_H = pms.Molecule(inp="data/heptane.gro")

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_density("output/dens_cyl_s.obj")
    sample.init_gyration("output/gyr_cyl_s.obj")
    sample.init_angle("output/angle_cyl.obj", [0, 3])
    sample.sample(is_parallel=False)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_density("output/dens_cyl_s_const_a.obj", bin_const_A=True)
    sample.sample(is_parallel=False)
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
    sample.init_diffusion_mc("output/diff_mc_cyl_s.obj", len_step=[1, 2, 5, 10, 20, 30, 40, 50, 100, 200, 250, 300, 350])
    sample.sample(is_parallel=False)

    sample = pa.Sample([6.00035, 6.00035, 19.09191], "data/traj_box.xtc", mol_H)
    sample.init_diffusion_mc("output/diff_mc_box.obj", len_step=[1, 2, 5, 10, 20, 30, 40, 50])
    sample.sample(shift=[0, 0, 3.3], is_parallel=False, is_pbc=True)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_density("output/dens_cyl_p.obj")
    sample.init_gyration("output/gyr_cyl_p.obj")
    sample.init_diffusion_bin("output/diff_cyl_p.obj")
    sample.sample(is_parallel=False, is_pbc=False)

    sample = pa.Sample("data/pore_system_cylinder_new.yml", "data/traj_cylinder.xtc", mol_B)
    sample.init_diffusion_mc("output/diff_mc_cyl_p.obj", len_step=[1, 2, 5, 10, 20, 30, 40])
    sample.sample(is_parallel=False, is_pbc=True)

    # Pre-compute a small MC output for file_to_text and other tests
    model = pa.CosineModel("output/diff_mc_cyl_s.obj", 6, 10)
    model._len_step = [10, 20, 40]
    pa.MC().run(model, "output/mc_fixture.yml", nmc_eq=500, nmc=500, is_print=False, is_parallel=False)

    return True
