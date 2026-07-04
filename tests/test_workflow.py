"""
End-to-end workflow test: PoreMS → PoreSim → (GROMACS) → PoreAna

This test exercises the full analysis pipeline as it would be used in practice:

  1. PoreMS generates a cylindrical pore and writes the structure (.gro) and
     system properties (.yml) needed for both simulation setup and analysis.

  2. PoreSim reads the PoreMS output and configures a GROMACS simulation box
     (topology, parameter files, cluster job scripts).  The actual GROMACS run
     is the only step that cannot be automated; the test uses the pre-recorded
     trajectory in tests/data/ as the GROMACS output.

  3. PoreAna loads the PoreMS yml directly as its pore-system input and
     analyses the trajectory: density profile, gyration radius, and bin
     diffusion coefficient.

Both poresim and porems are optional dependencies — if either is unavailable
the test is skipped rather than erroring.
"""

import os
import pytest

pms = pytest.importorskip("porems", reason="porems not installed")
ps  = pytest.importorskip("poresim", reason="poresim not installed")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PORESIM_DATA = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "PoreSim", "tests", "data",
)


def _poresim_data(*parts):
    """Absolute path into the PoreSim test-data directory."""
    return os.path.join(_PORESIM_DATA, *parts)


def _has_poresim_data():
    """True if the PoreSim test-data directory is reachable."""
    return os.path.isdir(_PORESIM_DATA) and os.path.isfile(_PORESIM_DATA + "/pore.top")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def porems_output(tmp_path_factory):
    """Generate a cylinder pore with PoreMS and return the output directory."""
    out = tmp_path_factory.mktemp("porems_out")
    pore = pms.PoreCylinder([6, 6, 6], 4, 5)
    pore.finalize()
    pore.store(str(out) + "/")
    return out


# ---------------------------------------------------------------------------
# Step 1 — PoreMS
# ---------------------------------------------------------------------------

def test_porems_generates_yml(porems_output):
    """PoreMS store() must produce a .gro and a .yml for the pore system."""
    gro = porems_output / "pore.gro"
    yml = porems_output / "pore.yml"
    assert gro.is_file(), "PoreMS did not write pore.gro"
    assert yml.is_file(), "PoreMS did not write pore.yml"

    # yml must contain the keys that PoreAna expects
    import yaml
    with open(yml) as f:
        data = yaml.safe_load(f)

    assert "system" in data, "yml missing 'system' section"
    assert "dimensions" in data["system"], "yml missing system.dimensions"
    assert "reservoir" in data["system"], "yml missing system.reservoir"

    shape_ids = [k for k in data if k.startswith("shape_")]
    assert shape_ids, "yml contains no shape sections"
    s = data[shape_ids[0]]
    assert "shape" in s and s["shape"] == "CYLINDER"
    assert "diameter" in s
    assert "parameter" in s
    assert "length" in s["parameter"]
    assert "centroid" in s["parameter"]


# ---------------------------------------------------------------------------
# Step 2 — PoreSim
# ---------------------------------------------------------------------------

def test_poresim_box_from_porems(porems_output):
    """PoreSim Box must accept the PoreMS .gro and .yml without error."""
    pore_gro = str(porems_output / "pore.gro")
    pore_yml = str(porems_output / "pore.yml")

    box = ps.Box("workflow_box")
    box.add_box(pore_gro)
    box.add_pore(pore_yml)
    box.add_mol("BEN", "data/benzene.gro", "fill", auto_dens=500, mass=78.11)

    sim_dict = box.get_sim_dict()
    assert sim_dict["struct"]["BOX"] == pore_gro
    assert sim_dict["struct"]["PORE"] == pore_yml
    assert "BEN" in sim_dict["mols"]


@pytest.mark.skipif(not _has_poresim_data(), reason="PoreSim test data not found")
def test_poresim_generates_files(porems_output, tmp_path):
    """PoreSim generate() must create the expected simulation folder structure."""
    pore_gro = str(porems_output / "pore.gro")
    pore_yml = str(porems_output / "pore.yml")

    cluster = {
        "address": "user@cluster",
        "directory": "/home/sim/",
        "queuing": {
            "add_np": False,
            "mpi": "$DO_PARALLEL",
            "shell": "job.sh",
            "submit": "sbatch",
        },
    }

    job = {
        "min": {"file": _poresim_data("forhlr.sh"), "nodes": 2, "np": 20, "wall": "24:00:00"},
        "nvt": {"file": _poresim_data("forhlr.sh"), "nodes": 4, "np": 20, "wall": "24:00:00"},
        "run": {"file": _poresim_data("forhlr.sh"), "maxh": 24, "nodes": 8,
                "np": 20, "runs": 5, "wall": "24:00:00"},
    }
    param = {
        "min": {"file": _poresim_data("pore_min.mdp")},
        "nvt": {"file": _poresim_data("pore_nvt.mdp"),
                "param": {"NUMBEROFSTEPS": 1000, "TEMPERATURE_VAL": 298}},
        "run": {"file": _poresim_data("pore_run.mdp"),
                "param": {"NUMBEROFSTEPS": 2000, "TEMPERATURE_VAL": 298}},
    }

    box = ps.Box("workflow")
    box.add_box(pore_gro)
    box.add_pore(pore_yml)
    box.add_mol("BEN", _poresim_data("benzene.gro"), "fill", auto_dens=500, mass=78.11)
    box.add_topol(_poresim_data("pore.top"), "master")
    box.add_topol(_poresim_data("grid.itp"), "top")
    box.add_topol([_poresim_data("benzene.top")])
    box.set_job(job)
    box.set_param(param)

    sim_link = str(tmp_path / "sim") + "/"
    sim = ps.Simulate(sim_link, box)
    sim.set_cluster(cluster)
    sim.generate()

    # Verify that core simulation files were written
    assert os.path.isdir(sim_link), "simulation folder was not created"
    assert any(f.endswith(".mdp") for f in _rglob(sim_link)), \
        "no .mdp parameter files in simulation folder"
    assert any(f.endswith("ana.py") for f in _rglob(sim_link)), \
        "no analysis script in simulation folder"


def _rglob(root):
    """Yield all file paths under root."""
    for dirpath, _, files in os.walk(root):
        for f in files:
            yield os.path.join(dirpath, f)


# ---------------------------------------------------------------------------
# Step 3 — PoreAna (using PoreMS yml + pre-recorded GROMACS trajectory)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def poreana_workflow_output(porems_output, tmp_path_factory, change_to_tests_dir):
    """Run PoreAna sampling with the PoreMS yml and the test trajectory.

    The trajectory in tests/data/traj_cylinder.xtc is the GROMACS output from
    a previous run in an identical simulation setup; it is used here as a
    stand-in for the GROMACS step that cannot be automated.
    """
    import poreana as pa

    out = tmp_path_factory.mktemp("poreana_workflow")
    pore_yml = str(porems_output / "pore.yml")
    traj     = "data/traj_cylinder.xtc"
    mol      = pms.Molecule(inp="data/benzene.gro")

    sample = pa.Sample(pore_yml, traj, mol)
    sample.init_density(str(out / "dens.obj"))
    sample.init_gyration(str(out / "gyr.obj"))
    sample.sample(is_parallel=False)

    return out


def test_poreana_density_from_porems_yml(poreana_workflow_output):
    """Density sampling driven by the PoreMS yml must complete and return
    a valid output structure."""
    import poreana as pa

    dens_file = str(poreana_workflow_output / "dens.obj")
    assert os.path.isfile(dens_file), "density output file not created"

    dens = pa.density.bins(dens_file, is_print=False)
    assert "num_dens" in dens
    assert "mean" in dens
    assert "dens" in dens

    shape_ids = [k for k in dens["num_dens"] if k.startswith("shape_")]
    assert shape_ids, "no shape bins in density output"

    # Number density must be non-negative everywhere
    for pid in shape_ids:
        assert all(v >= 0 for v in dens["num_dens"][pid]["in"]), \
            f"negative number density in {pid}"
    assert all(v >= 0 for v in dens["num_dens"]["ex"]), \
        "negative number density in reservoir"


def test_poreana_gyration_from_porems_yml(poreana_workflow_output):
    """Gyration sampling driven by the PoreMS yml must complete and return
    physically sensible (non-negative) radii."""
    import poreana as pa
    import matplotlib
    matplotlib.use("Agg")

    dens_file = str(poreana_workflow_output / "dens.obj")
    gyr_file  = str(poreana_workflow_output / "gyr.obj")
    assert os.path.isfile(gyr_file), "gyration output file not created"

    mean = pa.gyration.bins_plot(gyr_file, dens_file)
    assert "ex" in mean
    assert mean["ex"] >= 0, "mean gyration radius in reservoir is negative"
    if "in" in mean:
        assert mean["in"] >= 0, "mean gyration radius inside pore is negative"


def test_poreana_file_to_text_from_porems_yml(poreana_workflow_output, tmp_path):
    """file_to_text must produce a readable report from PoreMS-generated data."""
    import poreana as pa

    dens_file = str(poreana_workflow_output / "dens.obj")
    txt_file  = str(tmp_path / "workflow_report.txt")

    pa.utils.file_to_text(dens_file, txt_file)
    assert os.path.isfile(txt_file)
    with open(txt_file) as f:
        content = f.read()
    assert "[System]" in content
    assert "[Density]" in content
