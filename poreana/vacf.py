import math
import numpy as np
from scipy.integrate import cumulative_trapezoid

import poreana.utils as utils

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

    num_new_time_origins = sum([data[bin]["density"] for bin in range(bin_num)]) / num_res
    print("Averaged over ", num_new_time_origins, " new time origins.")

    avg_res_per_bin = [data[bin]["density"] / num_new_time_origins for bin in range(bin_num)]
    print("Average number of residues per bin:", avg_res_per_bin)
    avg_res_per_bin = np.where(avg_res_per_bin == 0, np.nan, avg_res_per_bin)

    vacf_data = np.zeros((bin_num, num_res, corr_steps, 3))
    for bin in range(bin_num):
        vacf_data[bin] = np.array([data[bin][res_id] for res_id in data[bin] if isinstance(res_id, int)])
        vacf_data[bin] *= num_res / data[bin]["density"]
    
    integrated = np.zeros(vacf_data.shape)

    for bin in range(bin_num):
        for res_id in range(num_res):
            for dim in range(3):
                integrated[bin, res_id, :, dim] = cumulative_trapezoid(
                    vacf_data[bin, res_id, :, dim], dx=len_frame * sample_step, initial=0)

    return integrated
"""
       N   \langle v_i,l(0) v_i(t) \rangle      num_res * sum_vacf                         vacf * num_res  
erg = -------------------------------------- = ---------------------------------------- = -----------------
       N_l                                      avg_res_per_bin * num_new_time_origins     sum_res_per_bin 
"""