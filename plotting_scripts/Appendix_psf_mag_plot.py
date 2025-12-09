import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from psfMC.analysis.plotting import _get_trace
import glob


# Updated function to load overlapping indices
def get_overlapping_indices(indices_1, indices_2):
    indices_1 = np.load(indices_1)
    print(indices_1)
    indices_2 = np.load(indices_2)
    print(indices_2)

    common_values, indices_in_1, indices_in_2 = np.intersect1d(indices_1, indices_2, return_indices=True)
    print(common_values)

    return indices_in_1, indices_in_2

def load_results(keys_to_load, folder, indices):
    # Load results from .npy files
    results = {}
    print("----------")
    print(folder)
    for key in keys_to_load:  # Assuming you have a list of keys to load
        # print("Loading ", f"/fred/oz183/sberger/paper_2_obs_bias/src/plotting_scripts/{folder}/{key}.npy")
        if "truth" in key:
            results[key] = np.load(f"/fred/oz183/sberger/paper_2_obs_bias/src/plotting_scripts/{folder}/{key}.npy", allow_pickle=True)
            print(key)
            print(results[key][indices])
    return results

# Example Usage
if __name__ == "__main__":
    mag_22psf = "22psf_collated_data"
    mag_17psf = "17psf_collated_data"
    mag_19p5psf = "19p5psf_collated_data"

    quasar_name = "J159"
    filter_type = "356W"
    keys_to_load = [
        f"observed_mag_{quasar_name}_{filter_type}",
        f"observed_mag_sigs_{quasar_name}_{filter_type}",
        f"truth_mag_wq_{quasar_name}_{filter_type}",
        f"truth_mag_sigs_wq_{quasar_name}_{filter_type}",
    ]

    indices_file_22 = mag_22psf + "/indices_used_J159_356W.npy"
    indices_file_19p5 = mag_19p5psf + "/indices_used_J159_356W.npy"

    indices_file_17 = mag_17psf + "/indices_used_J159_356W.npy"

    mag_19p5psf_indices, mag_17psf_indices = get_overlapping_indices(indices_file_19p5, indices_file_17)
    load_results(keys_to_load, mag_19p5psf, mag_19p5psf_indices)
    load_results(keys_to_load, mag_17psf, mag_17psf_indices)
    # print(mag_22psf_indices)