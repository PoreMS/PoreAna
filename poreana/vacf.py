import math
import numpy as np
from scipy.integrate import cumulative_trapezoid, cumulative_simpson

import poreana.utils as utils

def vacf(link_data):
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

    box = sample["box"]["length"]

    directions = [i for i in range(3) if i!= direction]
    surface = np.array(box[directions[0]]*box[directions[1]])
    distances = np.diff(bins, axis=0)
    volumes = surface * distances

    num_new_time_origins = sum([data[bin]["density"] for bin in range(bin_num)]) / num_res

    avg_res_per_bin = [data[bin]["density"] / num_new_time_origins for bin in range(bin_num)]

    # density = [avg_res_per_bin[bin] / volumes[bin] for bin in range(bin_num)]

    vacf_data = np.zeros((bin_num, num_res, corr_steps, 3))
    for bin in range(bin_num):
        vacf_data[bin] = np.array([data[bin][res_id] for res_id in data[bin] if isinstance(res_id, int)]) / num_new_time_origins
        vacf_data[bin] *= num_res / avg_res_per_bin[bin]

    integrated_1 = np.zeros(vacf_data.shape)
    integrated_2 = np.zeros(vacf_data.shape)

    for bin in range(bin_num):
        for res_id in range(num_res):
            for dim in range(3):
                integrated_1[bin, res_id, :, dim] = cumulative_simpson(vacf_data[bin, res_id, :, dim], dx=len_frame, initial=0)
                integrated_2[bin, res_id, :, dim] = cumulative_trapezoid(vacf_data[bin, res_id, :, dim], dx=len_frame, initial=0)