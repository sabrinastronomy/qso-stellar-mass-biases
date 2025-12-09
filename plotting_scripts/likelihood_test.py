"""
Extract psfmc output for one galaxy
"""
num_samples = 2000
from psfMC.analysis import plot_hist, corner_plot
from psfMC import model_galaxy_mcmc, load_database
from psfMC.analysis.plotting import _get_trace
from matplotlib import pyplot as plt
import os
import numpy as np
import sys
from scipy.stats import norm

from astropy.io import fits
# sys.path.insert(1, '/fred/oz183/sberger/paper_1_bluetides/main_scripts')
# from get_gal_lums import GetGalLums
from scipy.stats import gaussian_kde

direc = "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES"

truth_mag_top = 100
bic = False

### defaults for paper plots
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})

# load in all indices selected for sample
indices_1 = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/selected_indices/selected_indices_sample_-25_-20.npy")
indices_2 = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/selected_indices/selected_indices_sample_-21_-18.npy")
indices_3 = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/selected_indices/all_indices_sample_-25_-24.npy")
indices_4 = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/selected_indices/FULL_selected_indices_sample_23.5_24.5.npy")
indices_5 = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/J159_photo_cropped_selected_indices_sample_23.48_24.48.npy")

indices = np.concatenate((indices_1, indices_2, indices_3, indices_5))
indices = indices.flatten()

# initialize arrays of data we want to save
residuals = []
observed_mag = []
observed_mag_wq_likelihoods = []
observed_mag_sigs = []
photo_observed_mag_arr = []
photo_truth_mag_arr = []
truth_mag_wq = []
truth_mag_wq_likelihoods = []
truth_mag_sigs_wq = []
bic_value = []
likelihoods_wo_quasar = []
indices_used = []

# quasar magnitude and other run variables
mag = 19.5
mag = str(np.round(mag, 1))
mag = mag.replace(".", "p")
filter_type = "356W"
ss = 2
exp_time = 10000

count = 0
for gal in indices:
    bic_file = f"{direc}/all_out_files_bic_{filter_type}_{mag}_full_quasar_fits_files_6p5_exp_{exp_time}_FOV_25.0_samp_{ss}/BIC_out_{gal}_db.fits"
    db_file = f"{direc}/all_out_files_{filter_type}_{mag}_full_quasar_fits_files_6p5_exp_{exp_time}_FOV_25.0_samp_{ss}/quasar_out_{gal}_db.fits"
    db_file_woq = f"{direc}/all_out_files_{filter_type}_wo_quasar_fits_files_6p5_exp_{exp_time}_FOV_25.0_samp_{ss}/without_quasar_out_{gal}_db.fits"

    # be careful with sorting with files below
    gal_only_file_psf_subtracted = f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_out_files_{filter_type}_{mag}_full_quasar_fits_files_6p5_exp_{exp_time}_FOV_25.0_samp_{ss}/quasar_out_{gal}_point_source_subtracted.fits"
    gal_only_file_no_quasar = f"/fred/oz183/sberger/paper_2_obs_bias/fit_files_gal_only/{filter_type}_wo_quasar_fits_files_6p5_exp_{exp_time}_FOV_25.0_samp_{ss}/exp_time_{exp_time}_without_quasar_ss_{ss}_not_smoothed_{gal}__with_background_JWST.NIRCAM.F{filter_type}.fits"
    try:
        db = load_database(db_file)
        if bic:
            db_bic = load_database(bic_file)
        db_woq = load_database(db_file_woq)

        if bic:
            hdul_bic = fits.open(bic_file)
        hdul_quasar = fits.open(db_file)
        hdul_wo_quasar = fits.open(db_file_woq)

        if bic:
            # Access the primary HDU (Header Data Unit)
            ln_bic = hdul_bic[1].data["lnprobability"]
            MCITER_bic = num_samples

            ln_quasar = hdul_quasar[1].data["lnprobability"]
            MCITER_quasar = num_samples

            ln_wo_quasar = hdul_wo_quasar[1].data["lnprobability"]
            likelihoods_wo_quasar.append(np.min(ln_wo_quasar))
            bic_wo_quasar = 1 * np.log(MCITER_bic) - 2 * np.min(ln_bic)
            bic_quasar = 8 * np.log(MCITER_quasar) - 2 * np.min(ln_quasar)

            diff = bic_wo_quasar - bic_quasar
        indices_used.append(gal)
        print("opened")
        # if diff > 10:
        # os.system(
        #         f"python image_check.py --quasar_mag=19p5 --filter_type=356W --ss=2 --exp_time=10000 --quasar_num={gal} --direc=debug_gt_100k")
    except:
        count += 1
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
    average_mag_galaxy_obs = np.average(trace)
    observed_mag.append(average_mag_galaxy_obs)
    observed_mag_wq_likelihoods.append(trace)
    observed_mag_sigs.append(np.std(trace))
    if bic:
        bic_value.append(diff)
    
    trace_name = "1_Sersic_mag"
    trace = _get_trace(trace_name, db_woq)
    average_mag_galaxy_true = np.average(trace)
    truth_mag_wq.append(average_mag_galaxy_true)
    truth_mag_wq_likelihoods.append(trace)
    truth_mag_sigs_wq.append(np.std(trace))
    if average_mag_galaxy_true < truth_mag_top:
        truth_mag_top = average_mag_galaxy_true
        obs_mag_top = average_mag_galaxy_obs
        name_of_gal = gal_only_file_no_quasar

np.save(f"observed_mag_wq_likelihoods_{mag}_{filter_type}.npy", observed_mag_wq_likelihoods) # photometric magnitudes
np.save(f"truth_mag_wq_likelihoods_{mag}_{filter_type}.npy", truth_mag_wq_likelihoods) # photometric magnitudes

truth_mag_wq = np.asarray(truth_mag_wq)
# sorting_arr = np.argsort(truth_mag_wq)
# truth_mag_wq = truth_mag_wq[sorting_arr]
truth_mag_wq_likelihoods = np.squeeze(np.asarray(truth_mag_wq_likelihoods))
observed_mag_wq_likelihoods = np.squeeze(np.asarray(observed_mag_wq_likelihoods))
# truth_mag_wq_likelihoods = truth_mag_wq_likelihoods[sorting_arr] # resorting by increasing mag so we can bin
# observed_mag_wq_likelihoods = observed_mag_wq_likelihoods[sorting_arr]
observed_mag = np.asarray(observed_mag)
# observed_mag = observed_mag[sorting_arr] # just mean of posterior
indices_used = np.asarray(indices_used)
# indices_used = indices_used[sorting_arr]
for i in range(len(truth_mag_wq)):
    plt.hist(observed_mag_wq_likelihoods[i], label=truth_mag_wq[i], density=True)
plt.legend()
plt.xlabel("Observed Magnitudes")
plt.ylabel("Probability")
plt.savefig("../image_pngs/bins.png")
plt.close()


np.save(f"likelihood_test_PHOTO_observed_mag_SAME_MAG_{mag}_{filter_type}.npy", photo_observed_mag_arr) # photometric magnitudes
np.save(f"likelihood_test_PHOTO_truth_mag_SAME_MAG_{mag}_{filter_type}.npy", photo_truth_mag_arr) # photometric magnitudes
np.save(f"likelihood_test_observed_mag_SAME_MAG_{mag}_{filter_type}.npy", observed_mag)
np.save(f"likelihood_test_observed_mag_sigs_SAME_MAG_{mag}_{filter_type}.npy", observed_mag_sigs)
np.save(f"WOQ_likelihood_test_observed_mag_SAME_MAG_{mag}_{filter_type}.npy", truth_mag_wq)
np.save(f"WOQ_likelihood_test_observed_mag_sigs_SAME_MAG_{mag}_{filter_type}.npy", truth_mag_sigs_wq)
np.save(f"indices_used_{mag}_{filter_type}.npy", indices_used)

mag_arr = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/full_6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()


# main observed vs truth galaxy magnitude plots
fig, ax1 = plt.subplots()
#
# bic_values_mask = np.asarray(bic_value) < 10
# bic_values_mask_good = np.invert(bic_values_mask)

sc = ax1.scatter(truth_mag_wq, observed_mag, c=np.asarray(truth_mag_wq) - np.asarray(observed_mag))
sc = ax1.scatter(truth_mag_wq, observed_mag, c="red", marker='x')
sc = ax1.scatter(truth_mag_wq, observed_mag,c=np.asarray(truth_mag_wq) - np.asarray(observed_mag))

ax1.set_xlabel(rf'$\rm m_{gal, truth}$ (F356W, 10ks)')
ax1.set_ylabel(rf'$\rm m_{gal, obs}$ (F356W, 10ks)')  # Label for the bottom x-axis

# Add a colorbar
cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
cbar.set_label(rf'Residuals $(\rm m_{gal, truth} - m_{gal, obs})$')

# Add a horizontal line labeled as "EIGER quasar J159–02"
y_value = 23.98  # Horizontal line position
errors = 0.16  # Size of the shaded region around the line

# # Add a horizontal line
# ax1.axhline(y=y_value, color='red', linestyle='--', label='EIGER quasar J159–02')
#
# # Add shading around the horizontal line with a size of 0.16 on either side
# ax1.axhspan(y_value - errors, y_value + errors, color='red', alpha=0.3)
#
# # Annotate the horizontal line with the label
# ax1.text(ax1.get_xlim()[0] * 1.1, y_value - 0.1, 'EIGER quasar J159–02',
#          color='red', verticalalignment='top')

# Reverse the x-axis and y-axis direction
ax1.invert_xaxis()
ax1.invert_yaxis()

# Or you can set it to a different range if needed
# ax2.set_xlim([min_value, max_value])
plt.savefig("../image_pngs/residuals.png", dpi=500)
plt.close()
exit()

# Main plot
fig, ax1 = plt.subplots()
print(indices_used)
sc = ax1.scatter(mag_arr[indices_used], observed_mag, c=mag_arr[indices_used] - np.asarray(observed_mag))
# sc = ax1.scatter(mag_arr[bic_values_mask], observed_mag[bic_values_mask], c="red", marker='x')

ax1.set_xlabel('photometric truth gal magnitude (F356W)')
ax1.set_ylabel('observed gal magnitude (F356W, 10ks)')  # Label for the bottom x-axis

# Add a colorbar
cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
cbar.set_label('residuals (photo truth gal mag - observed gal mag)')

# Add a horizontal line labeled as "EIGER quasar J159–02"
y_value = 23.98  # Horizontal line position
errors = 0.16  # Size of the shaded region around the line

# Add a horizontal line
ax1.axhline(y=y_value, color='red', linestyle='--', label='EIGER quasar J159–02')

# Add shading around the horizontal line with a size of 0.16 on either side
ax1.axhspan(y_value - errors, y_value + errors, color='red', alpha=0.3)

# Annotate the horizontal line with the label
ax1.text(ax1.get_xlim()[0] * 1.1, y_value - 0.1, 'EIGER quasar J159–02',
         color='red', verticalalignment='top')

# Reverse the x-axis and y-axis direction
ax1.invert_xaxis()
ax1.invert_yaxis()

# Or you can set it to a different range if needed
# ax2.set_xlim([min_value, max_value])
plt.savefig("../image_pngs/residuals_photo.png", dpi=500)
plt.close()

# Main plot
fig, ax1 = plt.subplots()
sc = ax1.scatter(mag_arr[indices_used], observed_mag, c=mag_arr[indices_used] - np.asarray(observed_mag))
# sc = ax1.scatter(mag_arr[bic_values_mask], observed_mag[bic_values_mask], c="red", marker='x')

ax1.set_xlabel('photometric truth gal magnitude (F356W)')
ax1.set_ylabel('observed gal magnitude (F356W, 10ks)')  # Label for the bottom x-axis

# Add a colorbar
cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
cbar.set_label('residuals (photo truth gal mag - observed gal mag)')

# Add a horizontal line labeled as "EIGER quasar J159–02"
y_value = 23.98  # Horizontal line position
errors = 0.16  # Size of the shaded region around the line

# Add a horizontal line
ax1.axhline(y=y_value, color='red', linestyle='--', label='EIGER quasar J159–02')

# Add shading around the horizontal line with a size of 0.16 on either side
ax1.axhspan(y_value - errors, y_value + errors, color='red', alpha=0.3)

# Annotate the horizontal line with the label
ax1.text(ax1.get_xlim()[0] * 1.1, y_value - 0.1, 'EIGER quasar J159–02',
         color='red', verticalalignment='top')

# Reverse the x-axis and y-axis direction
ax1.invert_xaxis()
ax1.invert_yaxis()

# Or you can set it to a different range if needed
# ax2.set_xlim([min_value, max_value])
plt.savefig("../image_pngs/residuals_photo.png", dpi=500)
plt.close()

fig, ax1 = plt.subplots()

sc = ax1.scatter(mag_arr[indices_used], truth_mag_wq, c=mag_arr[indices_used] - np.asarray(truth_mag_wq))
# sc = ax1.scatter(mag_arr[bic_values_mask], observed_mag[bic_values_mask], c="red", marker='x')
# Add a colorbar
cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
cbar.set_label('residuals (photo truth gal mag - observed gal mag)')

# Add a horizontal line labeled as "EIGER quasar J159–02"
y_value = 23.98  # Horizontal line position
errors = 0.16  # Size of the shaded region around the line

# Add a horizontal line
ax1.axhline(y=y_value, color='red', linestyle='--', label='EIGER quasar J159–02')

# Add shading around the horizontal line with a size of 0.16 on either side
ax1.axhspan(y_value - errors, y_value + errors, color='red', alpha=0.3)

# Annotate the horizontal line with the label
ax1.text(ax1.get_xlim()[0] * 1.1, y_value - 0.1, 'EIGER quasar J159–02',
         color='red', verticalalignment='top')

# Reverse the x-axis and y-axis direction
ax1.invert_xaxis()
ax1.invert_yaxis()
ax1.set_xlabel('photometric truth gal magnitude (F356W)')
ax1.set_ylabel('truth gal magnitude (F356W, 10ks)')  # Label for the bottom x-axis
plt.savefig("../image_pngs/sersic_vs_photo_shortened.png", dpi=500)
plt.close()


fig, ax1 = plt.subplots()
sc = ax1.scatter(np.arange(len(truth_mag_wq)), truth_mag_wq, label="Sersic mag")
sc = ax1.scatter(np.arange(len(truth_mag_wq)), mag_arr[indices_used], label="Photo mag")
plt.legend()
ax1.set_ylabel('gal magnitude (F356W)')
ax1.set_xlabel('gal num')  # Label for the bottom x-axis
plt.savefig("../image_pngs/mag_checks_shortened.png", dpi=500)
plt.close()
