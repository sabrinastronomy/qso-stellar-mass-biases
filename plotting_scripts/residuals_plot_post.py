import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import gaussian_kde


### defaults for paper plots
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})

# quasar magnitude and other run variables
mag = 19.5
mag = str(np.round(mag, 1))
mag = mag.replace(".", "p")
filter_type = "356W"
ss = 2
exp_time = 10000

photo_observed_mag_arr = np.load(f"likelihood_test_PHOTO_observed_mag_SAME_MAG_{mag}_{filter_type}_{exp_time}.npy")
photo_truth_mag_arr = np.load(f"likelihood_test_PHOTO_truth_mag_SAME_MAG_{mag}_{filter_type}_{exp_time}.npy")
observed_mag = np.load(f"likelihood_test_observed_mag_SAME_MAG_{mag}_{filter_type}_{exp_time}.npy")
observed_mag_sigs = np.load(f"likelihood_test_observed_mag_sigs_SAME_MAG_{mag}_{filter_type}_{exp_time}.npy")
truth_mag_wq = np.load(f"WOQ_likelihood_test_observed_mag_SAME_MAG_{mag}_{filter_type}_{exp_time}.npy")
truth_mag_sigs_wq = np.load(f"WOQ_likelihood_test_observed_mag_sigs_SAME_MAG_{mag}_{filter_type}_{exp_time}.npy")
indices_used = np.load(f"indices_used_{mag}_{filter_type}_{exp_time}.npy")

mag_arr = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/full_6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()


# main observed vs truth galaxy magnitude plots
fig, ax1 = plt.subplots()
#
# bic_values_mask = np.asarray(bic_value) < 10
# bic_values_mask_good = np.invert(bic_values_mask)

sc = ax1.scatter(truth_mag_wq, observed_mag, c=np.asarray(truth_mag_wq) - np.asarray(observed_mag))
sc = ax1.scatter(truth_mag_wq, observed_mag,c=np.asarray(truth_mag_wq) - np.asarray(observed_mag))

ax1.set_xlabel(r'$\mathrm{m_{\mathrm{gal, \ truth}}}$ (F356W, 10 ks)')
ax1.set_ylabel(r'$\mathrm{m_{\mathrm{gal, \ obs}}}$ (F356W, 10 ks)')

# ax1.set_xlabel(rf'$\rm m_{gal, truth}$ (F356W, 10ks)')
# ax1.set_ylabel(rf'$\rm m_{gal, obs}$ (F356W, 10ks)')  # Label for the bottom x-axis

# Add a colorbar
cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
cbar.set_label(r'$\mathrm{Residuals \ (m_{\mathrm{gal, \ truth}} - m_{\mathrm{gal, \ obs}})}$')

# Add a horizontal line labeled as "EIGER quasar J159–02"
y_value = 23.98  # Horizontal line position
errors = 0.16  # Size of the shaded region around the line

# # # Add a horizontal line
ax1.axhline(y=y_value, color='red', linestyle='--', label='EIGER quasar J159–02')
#
# # # Add shading around the horizontal line with a size of 0.16 on either side
ax1.axhspan(y_value - errors, y_value + errors, color='red', alpha=0.3)

# Define the red line position (approximate value from the plot)
red_line_obs_value = y_value
tolerance = errors  # Allowable range around the red line to capture nearby points

# Filter truth values near the red line
filtered_truth_values = truth_mag_wq[(observed_mag >= red_line_obs_value - tolerance) &
                             (observed_mag <= red_line_obs_value + tolerance)]

# Plot the histogram of filtered truth values
# Compute KDE of filtered values
kde = gaussian_kde(filtered_truth_values)
x_range = np.linspace(min(filtered_truth_values) -3, max(filtered_truth_values) + 3, 1000)
density = kde(x_range)
plt.plot(x_range, density + red_line_obs_value, color='red', linestyle='-', linewidth=2)  # Shift KDE up to match red line

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

