"""
Extract psfmc output for one galaxy
"""

from psfMC.analysis import plot_hist, corner_plot
from psfMC import model_galaxy_mcmc, load_database
from psfMC.analysis.plotting import _get_trace
from matplotlib import pyplot as plt
import os
import numpy as np
import sys
from astropy.io import fits
sys.path.insert(1, '/fred/oz183/sberger/paper_1_bluetides/main_scripts')
from get_gal_lums import GetGalLums

magnitudes_range = np.arange(17, 25, 0.1)
filter_types = ["150W", "356W"]
super_samples = [1, 2]

np.save(f"observed_mag_range_SAME_GAL_0.npy", magnitudes_range)

indices = [0, 1, 2, 3, 4, 5]
residuals = []
observed_mag = []
observed_mag_sigs = []

gal_lum_class = GetGalLums(filt=None, width_by_rad=False) # filter doesn't matter since not running full gal lum computation

for filter_type, ss in zip(filter_types, super_samples):
    for gal in indices:
        residuals = []
        observed_mag = []
        observed_mag_sigs = []
        photo_observed_mag_arr = []
        for mag in magnitudes_range:
            mag = str(np.round(mag, 1))
            mag = mag.replace(".", "p")
            db_file = f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_out_files_{filter_type}_{mag}_full_quasar_fits_files_6p5_exp_3100_FOV_25.0_samp_{ss}/quasar_out_{gal}_db.fits"
            gal_only_file_psf_subtracted = f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_out_files_{filter_type}_{mag}_full_quasar_fits_files_6p5_exp_3100_FOV_25.0_samp_{ss}/quasar_out_{gal}_point_source_subtracted.fits"
            gal_only_file_no_quasar = f"/fred/oz183/sberger/paper_2_obs_bias/fit_files_gal_only/{filter_type}_wo_quasar_fits_files_6p5_exp_3100_FOV_25.0_samp_{ss}/exp_time_3100_without_quasar_ss_{ss}_not_smoothed_{gal}__with_background_JWST.NIRCAM.F{filter_type}.fits"

            try:
                db = load_database(db_file)
                image_psf = fits.getdata(gal_only_file_psf_subtracted, ext=0)
                image_wo = fits.getdata(gal_only_file_no_quasar, ext=0)
            except:
                print("-----------------")
                print("gal_only_file_no_quasar", gal_only_file_no_quasar)
                print("gal_only_file_psf_subtracted", gal_only_file_psf_subtracted)
                print("db file", db_file)
                print("This didn't work:")
                print(f"mag {mag}")
                print(f"filter_type {filter_type}")
                print(f"ss {ss}")
                print(f"gal {gal}")
                print("-----------------")
                continue

            trace_name = "2_Sersic_mag"
            trace = _get_trace(trace_name, db)
            average_mag_galaxy = np.average(trace)
            observed_mag.append(average_mag_galaxy)
            observed_mag_sigs.append(np.std(trace))

            np.save(f"observed_mag_SAME_GAL_{gal}_{filter_type}.npy", observed_mag)
            np.save(f"observed_mag_sigs_SAME_GAL_{gal}_{filter_type}.npy", observed_mag_sigs)

            photo_observed_mag = gal_lum_class.FROM_IMAGE_get_apparent_magnitude(image_psf, gal)
            photo_truth_mag = gal_lum_class.FROM_IMAGE_get_apparent_magnitude(image_wo, gal)

            photo_observed_mag_arr.append(photo_observed_mag)
            np.save(f"PHOTO_observed_mag_SAME_GAL_{gal}_{filter_type}.npy", photo_observed_mag_arr) # photometric magnitudes
            np.save(f"PHOTO_truth_mag_SAME_GAL_{gal}_{filter_type}.npy", [photo_truth_mag]) # photometric magnitudes
