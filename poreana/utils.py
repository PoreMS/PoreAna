################################################################################
# Utils                                                                        #
#                                                                              #
"""Here popular basic methods are noted."""
################################################################################


import os
import time
import yaml
import pickle
import datetime
import numpy as np
import pandas as pd
import poreana as pa


def mkdirp(directory):
    """Create directory if it does not exist.

    Parameters
    ----------
    directory : string
        Directory name
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def column(data):
    """Convert given row list matrix into column list matrix.

    Parameters
    ----------
    data : list
        Row data matrix

    Returns
    -------
    data_col : list
        Column data matrix
    """
    return [list(col) for col in zip(*data)]


def tic():
    """MATLAB tic replica - return current time.

    Returns
    -------
    time : float
        Current time in seconds
    """
    return time.time()


def toc(t, message="", is_print=True):
    """MATLAB toc replica - return time difference to tic and alternatively
    print a message.

    Parameters
    ----------
    t : float
        Starting time - given by tic function
    message : string, optional
        Custom output message
    is_print : bool, optional
        True for printing an output message

    Returns
    -------
    time : float
        Time difference
    """
    if message:
        message += " - runtime = "

    t_diff = time.time() - t

    if is_print:
        print(message + "%6.3f" % t_diff + " s")

    return t_diff


def _to_python(obj):
    """Recursively convert numpy types to native Python for YAML serialisation."""
    if isinstance(obj, dict):
        return {k: _to_python(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [_to_python(v) for v in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def save(obj, link):
    """Save an object file or yaml file using pickle.

    Parameters
    ----------
    obj : Object
        Object to be saved
    link : string
        Specific link to save object (.obj or .yml)
    """
    ext = link.split(".")[-1]

    if ext == "yml":
        with open(link, "w") as f:
            f.write(yaml.dump(_to_python(obj)))
    elif ext == "obj":
        with open(link, "wb") as f:
            pickle.dump(obj, f)
    else:
        print("Wrong data type — supported formats: .obj, .yml")


def load(link, file_type=""):
    """Load pickled object files or yaml files.

    Parameters
    ----------
    link : string
        Specific link to load object
    file_type : string, optional
        Specify filetype - **obj** or **yml** — leave empty for automatic
        determination from the file extension

    Returns
    -------
    obj : Object
        Loaded object
    """
    ext = file_type if file_type else link.split(".")[-1]

    if ext == "yml":
        with open(link, "r") as f:
            return yaml.load(f, Loader=yaml.SafeLoader)
    elif ext == "obj":
        with open(link, "rb") as f:
            return pickle.load(f)
    else:
        print("Wrong data type — supported formats: .obj, .yml")
        return None


def file_to_text(link, link_output, link_dens=[]):
    """Convert an output data file to a human-readable text file.
    For bin-diffusion and gyration output, a density sampling file is
    additionally required via ``link_dens``.

    Parameters
    ----------
    link : string
        Link to output obj or yml file
    link_output : string
        Link to the output text file
    link_dens : string, optional
        Link to density output obj or yml file (required for diff_bin and gyr_bin)
    """
    data = load(link)

    ###############################
    # Further properties Function #
    ###############################
    if data["type"] == "gyr_bin":
        if not link_dens:
            print("Gyration calculation needs a density sampling file. Check documentation and set a link_dens.")
            return

        if "pore" in data:
            pore = data["pore"]
            first_shape = next(k for k in pore if k.startswith("shape"))
            data_system = [
                [["%.2f" % i for i in pore["box"]["dimensions"]]],
                [float(pore[first_shape]["diam"])],
                [float(pore["box"]["res"])],
                [pore[first_shape]["type"]],
            ]
            df_system = pd.DataFrame(data_system,
                                     index=["Box dimension (nm)", "Pore diameter (nm)", "reservoir (nm)", "type"],
                                     columns=["Value"])
        elif "box" in data:
            box = data["box"]["length"]
            df_system = pd.DataFrame([[["%.2f" % i for i in box]]],
                                     index=["Box dimension (nm)"], columns=["Value"])
        df_system = df_system.rename_axis("# Identifier", axis=1)

        inp_table = [[data["inp"]["bin_num"]], [data["inp"]["entry"]],
                     [data["inp"]["num_frame"]], [data["inp"]["mass"]]]
        df_inputs = pd.DataFrame(inp_table,
                                 index=["Bin number", "Entry", "Frame number", "Mass"],
                                 columns=["Value"])
        df_inputs = df_inputs.rename_axis("# Identifier", axis=1)

        gyr_in = pa.gyration.bins_plot(link, link_dens, intent="ex")
        gyr_pd = pd.DataFrame(gyr_in, index=["gyration"])
        gyr_pd = gyr_pd.rename_axis("# Identifier", axis=1)

        with open(link_output, "w") as f:
            f.write("# This file was created " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") + "\n")
            f.write("# Analysis created using PoreAna\n")
            f.write("# Object file : " + os.path.dirname(os.path.abspath(__file__)) + link + "\n\n")
            f.write("[System]\n")
            f.write(df_system.to_string())
            f.write("\n\n[Input]\n")
            f.write(df_inputs.to_string())
            f.write("\n\n[Gyration]\n")
            f.write(gyr_pd.to_string())

    ##########################
    # Bin Diffusion Function #
    ##########################
    elif data["type"] == "diff_bin":
        if not link_dens:
            print("Bin diffusion needs a density sampling file. Check documentation and set a link_dens.")
            return

        bins = pa.diffusion.bins(link, is_norm=True)
        dens = pa.density.bins(link_dens)
        mean = pa.diffusion.mean(bins, dens)

        if "pore" in data:
            pore = data["pore"]
            first_shape = next(k for k in pore if k.startswith("shape"))
            data_system = [
                [["%.2f" % i for i in pore["box"]["dimensions"]]],
                [float(pore[first_shape]["diam"])],
                [float(pore["box"]["res"])],
                [pore[first_shape]["type"]],
            ]
            df_system = pd.DataFrame(data_system,
                                     index=["Box dimension (nm)", "Pore diameter (nm)", "reservoir (nm)", "type"],
                                     columns=["Value"])
            df_system = df_system.rename_axis("# Identifier", axis=1)

        inp_table = [[data["inp"]["bin_num"]], [data["inp"]["entry"]],
                     [data["inp"]["num_frame"]], [data["inp"]["mass"]], [bins["is_norm"]]]
        df_inputs = pd.DataFrame(inp_table,
                                 index=["Bin number", "Entry", "Frame number", "Mass", "is_norm"],
                                 columns=["Value"])
        df_inputs = df_inputs.rename_axis("# Identifier", axis=1)

        df_data_dict = {}
        for pid in bins["width"]:
            df_data_dict[f"# Width ({pid})"] = bins["width"][pid][1:]
            df_data_dict[f"D ({pid})"] = bins["diff"][pid]
        df_data = pd.DataFrame(df_data_dict)

        with open(link_output, "w") as f:
            f.write("# This file was created " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") + "\n")
            f.write("# Analysis created using PoreAna\n")
            f.write("# Object file : " + os.path.dirname(os.path.abspath(__file__)) + link + "\n")
            f.write("# Units\n# Diffusion D(10^-9 m^2s^-1)\n\n")
            f.write("[System]\n")
            f.write(df_system.to_string())
            f.write("\n\n[Inputs]\n")
            f.write(df_inputs.to_string())
            f.write("\n\n[Diffusion]\n")
            for pid, val in mean.items():
                f.write("%.2f " % val + "* 10^-9 m^2s^-1  (" + pid + ")\n")
            f.write("\n\n[Bin Diffusion Profiles]\n")
            f.write(df_data.to_string(index=False))

    ####################
    # Density Function #
    ####################
    elif data["type"] == "dens_bin":
        if "pore" in data:
            pore = data["pore"]
            first_shape = next(k for k in pore if k.startswith("shape"))
            data_system = [
                [["%.2f" % i for i in pore["box"]["dimensions"]]],
                [float(pore[first_shape]["diam"])],
                [float(pore["box"]["res"])],
                [pore[first_shape]["type"]],
            ]
            df_system = pd.DataFrame(data_system,
                                     index=["Box dimension (nm)", "Pore diameter (nm)", "reservoir (nm)", "type"],
                                     columns=["Value"])
        elif "box" in data:
            df_system = pd.DataFrame([[["%.2f" % i for i in data["box"]["length"]]]],
                                     index=["Box dimension (nm)"], columns=["Value"])
        df_system = df_system.rename_axis("# Identifier", axis=1)

        inp_table = [
            [data["inp"]["bin_num"]],
            [data["inp"]["entry"]],
            [data["inp"]["num_frame"]],
            [bool(data["inp"]["remove_pore_from_res"])],
            ["%.2f" % data["inp"]["mass"]],
        ]
        df_inputs = pd.DataFrame(inp_table,
                                 index=["Bin number", "Entry", "Frame number",
                                        "Remove pore from reservoir", "Mass"],
                                 columns=["Value"])
        df_inputs = df_inputs.rename_axis("# Identifier", axis=1)

        dens = pa.density.bins(link)
        if "pore" in data:
            ads = pa.adsorption.calculate(link)
            df_ads = pd.DataFrame(ads).rename_axis("# Identifier", axis=1)
            shape_ids = [k for k in data["data"] if k.startswith("shape")]
            mean_rows = [
                [dens["mean"][pid]["in"], dens["dens"][pid]["in"]]
                for pid in shape_ids
            ] + [[dens["mean"]["ex"], dens["dens"]["ex"]]]
            mean_idx = [f"Density inside {pid}" for pid in shape_ids] + ["Density outside pore"]
            df_mean = pd.DataFrame(
                mean_rows, index=mean_idx,
                columns=["Density (#/nm^3)", "Density (kg/m^3)"],
            )
        else:
            df_mean = pd.DataFrame([[dens["mean"]["ex"], dens["dens"]["ex"]]],
                                   index=["Density box"],
                                   columns=["Density (#/nm^3)", "Density (kg/m^3)"])
        df_mean = df_mean.rename_axis("# Identifier", axis=1)

        data_dens = {"# Ex width": dens["sample"]["data"]["ex_width"],
                     "Density (Ex)": dens["num_dens"]["ex"]}
        if "pore" in data:
            first_pid = next(k for k in data["data"] if k.startswith("shape"))
            data_dens["In width"] = dens["sample"]["data"][first_pid]["in_width"][1:]
            data_dens["Density (In)"] = dens["num_dens"][first_pid]["in"]
        df_data = pd.DataFrame(data_dens)

        with open(link_output, "w") as f:
            f.write("# This file was created " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") + "\n")
            f.write("# Analysis created using PoreAna\n")
            f.write("# Object file : " + os.path.dirname(os.path.abspath(__file__)) + link + "\n\n")
            f.write("[System]\n")
            f.write(df_system.to_string())
            f.write("\n\n[Inputs]\n")
            f.write(df_inputs.to_string())
            if "pore" in data:
                f.write("\n\n[Adsorption]\n")
                f.write(df_ads.to_string(na_rep=""))
            f.write("\n\n[Density]\n")
            f.write(df_mean.to_string())
            f.write("\n\n[Density Profiles]\n")
            f.write(df_data.to_string(index=False))

    ###############
    # MC Function #
    ###############
    elif data["type"] == "mc":
        t = data["model"]["len_frame"]
        if "pore" in data:
            pore = data["pore"]
            first_shape = next(k for k in pore if k.startswith("shape"))
            sys_data = [
                [["%.2f" % i for i in pore["box"]["dimensions"]]],
                [float(pore[first_shape]["diam"])],
                [float(pore["box"]["res"])],
                [pore[first_shape]["type"]],
            ]
            df_system = pd.DataFrame(sys_data,
                                     index=["Box dimension (nm)", "Pore diameter (nm)", "reservoir (nm)", "type"],
                                     columns=["Value"])
        elif "box" in data:
            df_system = pd.DataFrame([[["%.2f" % i for i in data["box"]["length"]]]],
                                     index=["Box dimension (nm)"], columns=["Value"])
        df_system = df_system.rename_axis("# Identifier", axis=1)

        diff = pa.diffusion.mc_profile(link, is_plot=False, infty_profile=True)
        free_energy = pa.freeenergy.mc_profile(link, is_plot=False)

        profile_data = {"# Bins [nm]": diff[2]}
        for i in diff[1]:
            profile_data[f"D (t={t * i})"] = diff[1][i]
        profile_data["   D (t=∞)"] = diff[0]
        for i in free_energy[0].keys():
            profile_data["Free energy [-]"] = free_energy[0][i]
        df_data = pd.DataFrame(profile_data)

        df_inputs = pa.tables.mc_inputs(link, print_con=False).rename_axis("# Identifier", axis=1)
        df_model = pa.tables.mc_model(link, print_con=False).rename_axis("# Identifier", axis=1)
        df_results = pa.tables.mc_results(link, print_con=False).rename_axis("# Identifier", axis=1)

        with open(link_output, "w") as f:
            f.write("# This file was created " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S") + "\n")
            f.write("# Analysis created using PoreAna\n")
            f.write("# Object file : " + os.path.dirname(os.path.abspath(__file__)) + "/" + link + "\n")
            f.write("# Units\n# Diffusion D(10^-9 m^2s^-1)\n# Lag time  t(ps)\n\n")
            f.write("[System]\n")
            f.write(df_system.to_string())
            f.write("\n\n\n[Model Inputs]\n")
            f.write(df_model.to_string())
            f.write("\n\n\n[MC Inputs]\n")
            f.write(df_inputs.to_string())
            f.write("\n\n\n[MC Results]\n")
            f.write(df_results.to_string())
            f.write("\n\n[Profiles]\n")
            f.write(df_data.to_string(index=False))


def num_dens_to_mass_dens(dens):
    """Convert number density
    :math:`\\rho \\left(\\frac{\\text{#}}{\\text{nm}^3}\\right)`
    into mass density
    :math:`\\rho_{\\text{m}} \\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)`

    .. math::

        \\rho_{m}=\\frac{M \\cdot 10}{6.022} \\cdot \\rho.

    Parameters
    ----------
    dens : dictionary
        Dictionary returned from :func:`pa.density.bins`

    Returns
    -------
    mass_dens : dictionary
        Mass density over the simulation box
    """
    factor = dens["sample"]["inp"]["mass"] * 10 / 6.022
    return {
        "ex": [factor * d for d in dens["num_dens"]["ex"]],
        "in": [factor * d for d in dens["num_dens"]["in"]] if "pore" in dens["sample"] else [],
    }


def mumol_m2_to_mols(c, A):
    r"""Convert concentration in :math:`\frac{\mu\text{mol}}{\text{m}^2}`
    to number of molecules.

    .. math:: N = 0.6022 \cdot c \cdot A

    Parameters
    ----------
    c : float
        Concentration in :math:`\frac{\mu\text{mol}}{\text{m}^2}`
    A : float
        Surface in :math:`\text{nm}^2`

    Returns
    -------
    N : float
        Number of molecules
    """
    return 0.6022 * c * A


def mols_to_mumol_m2(N, A):
    r"""Convert number of molecules to concentration in
    :math:`\frac{\mu\text{mol}}{\text{m}^2}`.

    .. math:: c = \frac{N}{0.6022 \cdot A}

    Parameters
    ----------
    N : int
        Number of molecules
    A : float
        Surface in :math:`\text{nm}^2`

    Returns
    -------
    c : float
        Concentration in :math:`\frac{\mu\text{mol}}{\text{m}^2}`
    """
    return N / 0.6022 / A


def mmol_g_to_mumol_m2(c, SBET):
    r"""Convert concentration in :math:`\frac{\text{mmol}}{\text{g}}`
    to :math:`\frac{\mu\text{mol}}{\text{m}^2}`.

    .. math:: c_A = \frac{c_g}{S_\text{BET}} \cdot 10^3

    Parameters
    ----------
    c : float
        Concentration in :math:`\frac{\text{mmol}}{\text{g}}`
    SBET : float
        Material surface in :math:`\frac{\text{m}^2}{\text{g}}`

    Returns
    -------
    c : float
        Concentration in :math:`\frac{\mu\text{mol}}{\text{m}^2}`
    """
    return c / SBET * 1e3


def mmol_l_to_mols(c, V):
    r"""Convert concentration in :math:`\frac{\text{mmol}}{\text{l}}`
    to number of molecules.

    .. math:: N = 6.022 \cdot 10^{-4} \cdot c \cdot V

    Parameters
    ----------
    c : float
        Concentration in :math:`\frac{\text{mmol}}{\text{l}}`
    V : float
        Volume in :math:`\text{nm}^3`

    Returns
    -------
    N : float
        Number of molecules
    """
    return 6.022e-4 * c * V


def mols_to_mmol_l(N, V):
    r"""Convert number of molecules to concentration in
    :math:`\frac{\text{mmol}}{\text{l}}`.

    .. math:: c = \frac{N}{6.022 \times 10^{-4} \cdot V}

    Parameters
    ----------
    N : int
        Number of molecules
    V : float
        Volume in :math:`\text{nm}^3`

    Returns
    -------
    c : float
        Concentration in :math:`\frac{\text{mmol}}{\text{l}}`
    """
    return N / 6.022e-4 / V
