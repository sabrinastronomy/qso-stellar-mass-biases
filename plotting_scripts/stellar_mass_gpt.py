import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import quasar_dict

quasar_data = quasar_dict.quasar_data

# Configure matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})

# File paths
dataholder_file_name = "/fred/oz183/sberger/paper_1_bluetides/main_scripts/bluetides_glory/full_bluetides_6p5.pkl"
mags_115W = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_6_JWST.NIRCAM.F115W_mag_all_order_luminous_BH.npy").flatten()
mags_150W = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_6_JWST.NIRCAM.F150W_mag_all_order_luminous_BH.npy").flatten()
mags_200W = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_6_JWST.NIRCAM.F200W_mag_all_order_luminous_BH.npy").flatten()
mags_356W = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()
# indices = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/J1120_photo_cropped_selected_indices_sample_23.95_24.95.npy")
# indices_356W = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/J1120_photo_cropped_selected_indices_sample_23.95_24.95.npy")
# indices_150W = np.load("/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/J1120_photo_cropped_selected_indices_sample_23.95_24.95.npy")

# Stellar mass data
stellarMass = np.load('/fred/oz183/mmarshal//BlueTides/stellarMass_208.npy')[:108001]


# Colors for filters
filter_colors = {
    "115W": "#88CCEE",
    "150W": "#117733",
    "200W": "#FF9900",  # orange
    "356W": "#CC6677"
}

mag_dict = {
    "115W": mags_115W,
    "150W": mags_150W,
    "200W": mags_200W,
    "356W": mags_356W
}
zorder = 100
def plot_stellar_mass_with_quasars(mags_dict_curr, stellarMass, quasar_data, save_dir):
    """Plot stellar mass vs magnitude and overlay quasar data with vertical lines."""
    # plt.semilogy(mags_356W, stellarMass,
    #              marker='o', linestyle='', color="grey", label=f"all mag=356", alpha=0.1)


    # Overlay quasar galaxy magnitudes as vertical lines
    count = 0
    low_lim, up_lim = None, None
    # for quasar, data in quasar_data.items():
    for quasar in ["J2255+0251"]:
        data = quasar_data[quasar]
        for filt, mag_data in data["galaxy"].items():
            if mag_data is not None and filt in ["150W", "356W"]:
                mag, _ = mag_data

                mags_for_filter = mags_dict_curr[filt]
                mag_mask = (mags_for_filter > mag - 1) & (mags_for_filter < mag + 1)
                # print(len(mags_for_filter[mag_mask]))
                # Plot stellar mass vs. magnitude for each filter
                # if filt == "F356W":
                #     zorder = -100
                # else:
                #     zorder = 100
                plt.semilogy(mags_for_filter[mag_mask], stellarMass[mag_mask],
                             marker='o', linestyle='', color=filter_colors[filt], alpha=0.1, zorder=zorder)
                if count == 0:
                    low_lim, up_lim = plt.ylim()  # Get y-axis limits
                    count += 1
                plt.vlines(
                    mag, low_lim, up_lim,
                    colors=filter_colors[filt],
                    label=f"{filt}", alpha=1
                )
    # Plot settings
    plt.xlabel(r"$\rm m_{galaxy, JWST}~[mag]$")
    plt.ylabel(r"$\rm M_{*}~[M_{\odot}]$")
    plt.legend()
    plt.ylim(low_lim, up_lim)
    plt.xlim(24, 27)
    # Save plot
    plt.savefig(f"{save_dir}/stellar_mass_with_quasar_lines.png", dpi=150)

# Save directory
save_dir = "../posteriors"

# Plot with updated function
plot_stellar_mass_with_quasars(mag_dict, stellarMass, quasar_data, save_dir)
