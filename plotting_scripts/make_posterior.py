import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

### defaults for paper plots
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})

# Define the model function for fitting (e.g., exponential)
def model_func(x, a, b):
    return a * np.exp(b * x)

mag_arr = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()
mag_arr_no_rad = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()
obs_mag = np.load("likelihood_test_observed_mag_SAME_MAG_19p5_356W.npy")
truth_mag_wq = np.load("WOQ_likelihood_test_observed_mag_SAME_MAG_19p5_356W.npy")
obs_mag_posteriors = np.load(f"observed_mag_wq_likelihoods_19p5_356W.npy") # obs mag posteriors
print(len(obs_mag))
# correct indices
indices_used = np.load("indices_used_19p5_356W.npy")
num_bins = 10

# # Create 10 bins (equal-width)
bins = np.linspace(22, 25, num_bins + 1)
mag_truths_bin_centers = (bins[1:] + bins[:-1]) / 2
# Digitize the array elements into bins
bin_indices = np.digitize(truth_mag_wq, bins) - 1  # subtract 1 to make bins 0-indexed

# Create a list of indices for each bin
indices_per_bin = [np.where(bin_indices == i)[0].tolist() for i in range(num_bins)]
posterior_truth_mag_given_obs_mag = []
m_obs_eiger = 23.98
# mag_obs_trying = [m_obs_eiger, 23.6, 24.3, 24]
mag_obs_trying = [m_obs_eiger]

# Your data points
x_data = mag_truths_bin_centers
y_data = posterior_truth_mag_given_obs_mag



for m_obs in mag_obs_trying:
    posterior_truth_mag_given_obs_mag = []
    # get observed mag posteriors for each bin
    for j, indices_in_bin in enumerate(indices_per_bin):
        print("indices_in_bin")
        print(indices_in_bin)
        mag_obs_likelihood_posteriors_selected = obs_mag_posteriors[indices_in_bin].flatten()
        # plt.hist(mag_obs_likelihood_posteriors_selected, density=True)
        mu, std = norm.fit(mag_obs_likelihood_posteriors_selected)
        p = norm.pdf(m_obs, mu, std)
        posterior_truth_mag_given_obs_mag.append(p)

    # plt.savefig("../image_pngs/hist.png")
    # plt.close()
    posterior_truth_mag_given_obs_mag = np.where(np.isnan(posterior_truth_mag_given_obs_mag), 0, posterior_truth_mag_given_obs_mag)
    posterior_truth_mag_given_obs_mag = np.where(posterior_truth_mag_given_obs_mag < 1e-4, 0, posterior_truth_mag_given_obs_mag)
    mask = posterior_truth_mag_given_obs_mag != 0

    plt.title(f"m_obs = {m_obs}")
    plt.semilogy(mag_truths_bin_centers[mask], posterior_truth_mag_given_obs_mag[mask], 'o', c="black")
    plt.axvline(x=m_obs_eiger, color='black', linestyle='--', label='Observed EIGER Galaxy m = 23.98')
    print(mag_truths_bin_centers)
    print(posterior_truth_mag_given_obs_mag)

    # # Create the interpolation function (cubic for smoothness)
    interpolation_func = interp1d(mag_truths_bin_centers[mask], posterior_truth_mag_given_obs_mag[mask], kind='cubic', fill_value="extrapolate")

    # Generate new x values for a smooth curve
    x_smooth = np.linspace(min(mag_truths_bin_centers), max(mag_truths_bin_centers), 100)
    y_smooth = interpolation_func(x_smooth)
    print(x_smooth)
    print(y_smooth)
    plt.plot(x_smooth, y_smooth)
    plt.xlabel("m_truth")
    plt.ylabel("dP/dm_truth | m_obs")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../image_pngs_debug/posterior_{m_obs}.png", dpi=300)
    plt.close()
