import matplotlib.pyplot as plt
import numpy as np
from psfMC import model_galaxy_mcmc, load_database
from psfMC.analysis import plot_hist, corner_plot
from psfMC.analysis.plotting import _get_trace
from astropy.io import fits
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

galaxies = [0, 1, 2, 3, 4, 5]
for i in galaxies:
    for filter_type, ss in zip(["150W", "356W"], [1, 2]):
    # for filter_type, ss in zip(["150W"], [1]):

        print(i, filter_type, ss)
        db_file = f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_out_files_{filter_type}_wo_quasar_fits_files_6p5_exp_3100_FOV_25.0_samp_{ss}/without_quasar_out_{i}_db.fits"
        db  = load_database(db_file)

        trace_name = "1_Sersic_mag"
        trace = _get_trace(trace_name, db)
        TRUTH_MAG = np.average(trace)
        if filter_type == "150W":
            truth_mag_150w = np.round(TRUTH_MAG, 1)
        if filter_type == "356W":
            truth_mag_356w = np.round(TRUTH_MAG, 1)
        observed_gal_mag = np.load(f"observed_mag_SAME_GAL_{i}_{filter_type}.npy")
        observed_gal_mag_sig = np.load(f"observed_mag_sigs_SAME_GAL_{i}_{filter_type}.npy")
        photo_observed_gal_mag = np.load(f"PHOTO_observed_mag_SAME_GAL_{i}_{filter_type}.npy")
        truth_photo_gal_mag = np.load(f"PHOTO_truth_mag_SAME_GAL_{i}_{filter_type}.npy")[0]
        length = len(observed_gal_mag)
        quasar_mag = np.load(f"observed_mag_range_SAME_GAL_0.npy")[:length]

        plt.errorbar(quasar_mag, TRUTH_MAG - observed_gal_mag, yerr=observed_gal_mag_sig, label=fr"{filter_type}, #{i}", fmt="o")
        plt.errorbar(quasar_mag, truth_photo_gal_mag - photo_observed_gal_mag, yerr=0, label=fr"{filter_type}, #{i} (photo)", fmt="o", alpha=0.5)

    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), bbox_transform=plt.gca().transAxes)
    plt.xlabel(r"quasar mag")
    plt.ylabel(r"truth-observed magnitude")  # (> 0, galaxy stellar mass overestimated, < 0 galaxy stellar mass underestimated)

    obs_150w_file_without_quasar = f"/fred/oz183/sberger/paper_2_obs_bias/fit_files_gal_only/150W_wo_quasar_fits_files_6p5_exp_3100_FOV_25.0_samp_1/exp_time_3100_without_quasar_ss_1_not_smoothed_{i}__with_background_JWST.NIRCAM.F150W.fits"
    image_data_150 = fits.getdata(obs_150w_file_without_quasar, ext=0)

    obs_356w_file_without_quasar = f"/fred/oz183/sberger/paper_2_obs_bias/fit_files_gal_only/356W_wo_quasar_fits_files_6p5_exp_3100_FOV_25.0_samp_1/exp_time_3100_without_quasar_ss_1_not_smoothed_{i}__with_background_JWST.NIRCAM.F356W.fits"
    image_data_356 = fits.getdata(obs_356w_file_without_quasar, ext=0)

    center = np.shape(image_data_356)[0] // 2
    offset = 25

    image_data_356 = image_data_356[center - offset:center + offset, center - offset:center + offset]
    image_data_150 = image_data_150[center - offset:center + offset, center - offset:center + offset]

    axins1 = inset_axes(plt.gca(), width="30%", height="20%", loc="upper right")
    axins1.imshow(image_data_356)
    axins1.axis('off')


    axins1.text(0.5, 0.95, r"$m_{356w}$: " + str(truth_mag_356w), color='white', fontsize=5, ha='center', va='center',
                transform=axins1.transAxes)


    plt.savefig(f"../image_pngs/gal_comp_{i}.png", dpi=2000)
    plt.close()

