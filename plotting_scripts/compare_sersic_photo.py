import numpy as np
import matplotlib.pyplot as plt

mag_arr = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()
mag_arr_no_rad = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()
obs_mag = np.load("likelihood_test_observed_mag_SAME_MAG_19p5_356W.npy")

# correct indices
indices_used = np.load("indices_used_19p5_356W.npy")
mag_arr = mag_arr[indices_used]
mag_sersic = np.load("WOQ_likelihood_test_observed_mag_SAME_MAG_19p5_356W.npy")

sorted_indices = np.argsort(indices_used)

# sorting
mag_sersic = mag_sersic[sorted_indices]
mag_arr = mag_arr[sorted_indices]
mag_arr_no_rad = mag_arr_no_rad[sorted_indices]
obs_mag = obs_mag[sorted_indices]

fig, ax1 = plt.subplots()

sc = ax1.scatter(mag_sersic, mag_arr, alpha=0.5, label="radial cropped photo mag")
sc = ax1.scatter(mag_sersic, mag_arr_no_rad, alpha=0.1, label="cropped to 76ckpc square aperture photo mag")
sc = ax1.plot(mag_sersic, mag_sersic, label="1:1", color="black", linewidth=3)

# cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
# cbar.set_label('residuals (photo - sersicn mag)')
plt.legend()
ax1.set_ylabel('w/o quasar sersic gal magnitude (F356W)')
ax1.set_xlabel('photo mag (F356W')  # Label for the bottom x-axis
plt.savefig("../image_pngs_debug/photo_comp_sersic.png", dpi=500)
plt.close()

fig, ax1 = plt.subplots()

sc = ax1.scatter(np.arange(len(mag_sersic)), mag_sersic, label="Sersic mag", alpha=0.5)
sc = ax1.scatter(np.arange(len(mag_arr)), mag_arr, label="Photo mag (cropped to radius)", alpha=0.5)
sc = ax1.scatter(np.arange(len(mag_arr_no_rad)), mag_arr_no_rad, label="Photo mag (cropped to one square aperture)", alpha=0.5)

# cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
# cbar.set_label('residuals (photo - sersicn mag)')
plt.legend()
ax1.set_ylabel('gal magnitude (F356W)')
ax1.set_xlabel('gal num')  # Label for the bottom x-axis
plt.savefig("../image_pngs_debug/mag_checks_shortened.png", dpi=500)
plt.close()

# Main plot
fig, ax1 = plt.subplots()
sc = ax1.scatter(mag_arr, obs_mag, c=mag_arr - np.asarray(obs_mag))
# sc = ax1.scatter(mag_arr[bic_values_mask], observed_mag[bic_values_mask], c="red", marker='x')

ax1.set_xlabel('cropped by rad photometric truth gal magnitude (F356W)')
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
plt.savefig("../image_pngs_debug/residuals_photo.png", dpi=500)
plt.close()


fig, ax1 = plt.subplots()

sc = ax1.scatter(mag_sersic, obs_mag, c=np.asarray(mag_sersic) - np.asarray(obs_mag))

ax1.set_xlabel('truth gal magnitude (F356W, 10ks)')
ax1.set_ylabel('observed gal magnitude (F356W, 10ks)')  # Label for the bottom x-axis

# Add a colorbar
cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
cbar.set_label('residuals (truth gal mag - observed gal mag)')

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
plt.savefig("../image_pngs_debug/residuals.png", dpi=500)
plt.close()