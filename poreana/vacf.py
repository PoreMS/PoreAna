import math
import numpy as np
from scipy.integrate import cumulative_trapezoid

import poreana.utils as utils


def density(link_data):
    """
    Calculate the density from given VACF data.
    The density is returned as the average number of residues per bin.

    Parameters
    ----------
    link_data : str
        The path to the data file containing the VACF data.
    
    Returns
    -------
    avg_res_per_bin : np.ndarray
        The average number of residues per bin, with shape (bin_num,).
    """
    sample = utils.load(link_data)
    data = sample["data"]
    inp = sample["inp"]
    num_res = inp["num_res"]

    num_new_time_origins = np.sum(data["density"]) / num_res
    avg_res_per_bin = data["density"] / num_new_time_origins
    return avg_res_per_bin


def vacf(link_data):
    """
    Calculate the velocity autocorrelation function (VACF) from the given data.
    The VACF is calculated using the cumulative trapezoid rule for integration.

    Parameters
    ----------
    link_data : str
        The path to the data file containing the VACF data.
    Returns
    -------
    integrated : np.ndarray
        The integrated VACF data, with shape (bin_num, num_res, corr_steps, 3).
    """
    sample = utils.load(link_data)

    data = sample["data"]

    inp = sample["inp"]
    num_frame = inp["num_frame"]
    mass = inp["mass"]
    bins = inp["bins"]
    len_correration = inp["len_correration"]
    new_time_origin = inp["new_time_origin"]
    sample_step = inp["sample_step"]
    len_frame = inp["len_frame"]
    bin_num = inp["bin_num"]
    direction = inp["direction"]
    corr_steps = inp["corr_steps"]
    new_time_origin_steps = inp["new_time_origin_steps"]
    num_res = inp["num_res"]

    num_new_time_origins = np.sum(data["density"]) / num_res
    print("Averaged over ", num_new_time_origins, " new time origins.")

    vacf_data = data["vacf_data"].copy() * num_res / data["density"][:, np.newaxis, np.newaxis, np.newaxis]
    
    integrated = cumulative_trapezoid(
        vacf_data, dx=len_frame * sample_step, axis=1, initial=0)
    
    # To make it more readable
    integrated = integrated.swapaxes(1, 2)

    return integrated
"""
       N   \langle v_i,l(0) v_i(t) \rangle      num_res * sum_vacf                         vacf * num_res  
erg = -------------------------------------- = ---------------------------------------- = -----------------
       N_l                                      avg_res_per_bin * num_new_time_origins     sum_res_per_bin 
"""