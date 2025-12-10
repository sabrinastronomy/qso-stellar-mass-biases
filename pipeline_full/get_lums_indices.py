
import matplotlib.pyplot as plt
import numpy as np

# mag_arr = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/full_6_JWST.NIRCAM.F356W_mag_all_order_luminous_BH.npy").flatten()
# mag_arr_uv = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/full_6_FAKE.FAKE.1500_mag_all_order_luminous_BH.npy").flatten()
# maddie_uv = np.load("/fred/oz183/sberger/paper_1_bluetides/main_scripts/maddie_lum_files/lum_dust_6p5_1500.npy")

gal_mag_dict = {
    "J0148+0600": {"F115W": [23.48, (0.24, 0.24)], "F200W": [23.51, (0.15, 0.15)], "F356W": [22.61, (0.07, 0.07)]},
    "J159-02": {"F115W": [24.83, (0.06, 0.06)], "F200W": [24.82, (0.23, 0.23)], "F356W": [23.98, (0.16, 0.16)]},
    "J1120+0641": {"F115W": [25.94, (0.20, 0.20)], "F200W": [25.48, (0.37, 0.37)], "F356W": [25.72, (0.47, 0.47)]},
    "J2236+0032": {"F150W": [25.12, (0.29, 0.29)], "F356W": [23.12, (0.20, 0.20)]},
    # "J2255": {"F356W": [24.58, (0.30, 0.30)]}
}


# gal_mag_dict = {
#     "J0844-0132": {"F356W": [25.70, (0.34, 0.34)], "F150W": [26.79, (None, None)]},
#     "J0844-0052": {"F356W": [25.57, (0.51, 0.51)], "F150W": [25.34, (0.25, 0.25)]},
#     "J0911+0152": {"F356W": [26.56, (0.30, 0.30)], "F150W": [27.61, (None, None)]},
#     "J0918+0139": {"F356W": [25.28, (0.25, 0.25)], "F150W": [26.03, (None, None)]},
#     "J1425-0015": {"F356W": [24.12, (0.29, 0.29)], "F150W": [25.78, (0.29, 0.29)]},
#     "J1512+4422": {"F356W": [23.10, (0.28, 0.28)], "F150W": [24.07, (0.17, 0.17)]},
#     "J1525+4303": {"F356W": [24.61, (0.24, 0.24)], "F150W": [25.60, (0.21, 0.21)]},
#     "J1146+0124": {"F356W": [24.15, (0.16, 0.16)], "F150W": [25.62, (0.20, 0.20)]},
#     "J1146-0005": {"F356W": [26.38, (None, None)], "F150W": [28.04, (None, None)]},
#     "J0217-0208": {"F356W": [23.69, (0.18, 0.18)], "F150W": [24.29, (0.13, 0.13)]},
#     "J2255+0251": {"F356W": [24.58, (0.30, 0.30)], "F150W": [26.3, (None, None)]}
# }


for quasar in gal_mag_dict:
    quasar_vals = gal_mag_dict[quasar].keys()
    for filter_type in quasar_vals:
        mag_arr_crop = np.load(
            f"/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_6_JWST.NIRCAM.{filter_type}_mag_all_order_luminous_BH.npy").flatten()
        mag_arr = mag_arr_crop.copy()
        mag = gal_mag_dict[quasar][filter_type][0]
        top_mag = np.round(mag + 2, 2)
        bottom_mag = np.round(mag - 2, 2)

        # Create a mask for values between top mag and bottom mag
        mask = (mag_arr > bottom_mag) & (mag_arr < top_mag)

        # Get the indices where the mask is True
        indices = np.where(mask)[0]
        # np.save(f"BY_WIDTH_{quasar}_{filter_type}_all_indices_sample_{bottom_mag}_{top_mag}.npy", indices)
        # print("TOTAL NUMBER OF POTENTIAL GALXIES")
        # print(len(indices))

        # Define the number of bins
        num_bins = 20  # Example: Create the num bins between -22 and -20
        selected_num_bins = 10

        if len(indices) < num_bins * selected_num_bins:
            selected_indices = indices
        else:
            bins = np.linspace(bottom_mag, top_mag, num_bins + 1)
            bin_indices = np.digitize(mag_arr[indices], bins)

            # Clip bin indices to be within valid range (1 to num_bins)
            bin_indices = np.clip(bin_indices, 1, num_bins)

            selected_indices = []

            for bin_num in range(1, num_bins + 1):
                # Find all indices that fall into this bin
                bin_mask = bin_indices == bin_num
                bin_members = indices[bin_mask]

                if len(bin_members) == 0:
                    continue  # Skip empty bins

                # Sample either all or up to selected_num_bins
                sample_size = min(len(bin_members), selected_num_bins)
                selected = np.random.choice(bin_members, size=sample_size, replace=False)

                selected_indices.extend(selected)

            # Optional: convert to array and deduplicate (should be unnecessary but safe)
            selected_indices = np.unique(selected_indices)


        np.save(f"indices_to_fit_direc/BY_WIDTH_{quasar}_{filter_type}_photo_cropped_selected_indices_sample_{bottom_mag}_{top_mag}.npy", selected_indices)

        print("===========")
        print(quasar)
        print(filter_type)
        print(len(selected_indices))
        print("===========")