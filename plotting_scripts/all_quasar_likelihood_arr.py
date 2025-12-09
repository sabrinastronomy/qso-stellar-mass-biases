import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from psfMC import model_galaxy_mcmc, load_database
from psfMC.analysis import plot_hist, corner_plot
from psfMC.analysis.plotting import _get_trace
import glob
import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 12})

class QuasarParser:
    def __init__(self, quasar_name, quasar_mag_dict, filter_type="356W", ss=2, exp_time=10000, num_samples=2000, psf_mag=17, photometric_only=False, bic=False):
        self.quasar_name_long = quasar_name
        self.quasar_mag_dict = quasar_mag_dict
        self.filter_type = filter_type
        self.ss = ss
        self.exp_time = exp_time
        self.num_samples = num_samples
        self.psf_mag = psf_mag
        self.bic = bic
        mag = self.quasar_mag_dict[self.quasar_name_long]["galaxy"][self.filter_type][0]
        try:
            top_mag = np.round(mag + 2, 2)
            bottom_mag = np.round(mag - 2, 2)
        except:
            top_mag = bottom_mag = None # test None

        pattern = f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/indices_to_fit_direc/" \
                  f"BY_WIDTH_{self.quasar_name_long}_{self.filter_type}_photo_cropped_selected_indices_sample_*_*.npy"

        files = glob.glob(pattern)
        if len(files) != 1:
            raise ValueError(f"Expected exactly one match, got {len(files)}: {files}")

        self.indices_full_length = len(np.load(files[0]))
        key_prefix = f"{self.quasar_name_long}_{self.filter_type}_{self.exp_time}"

        self.results = {
            f"observed_mag_{key_prefix}": [],
            f"observed_mag_likelihoods_{key_prefix}": [],
            f"observed_mag_sigs_{key_prefix}": [],
            f"truth_mag_wq_{key_prefix}": [],
            f"truth_mag_wq_likelihoods_{key_prefix}": [],
            f"truth_mag_sigs_wq_{key_prefix}": [],
            f"indices_used_{key_prefix}": []
        }

        if self.bic:
            self.results[f"BIC_{key_prefix}"] = []
        self.direc = "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES"
        self.photometric_only = photometric_only            
        
        self.process_quasar()
        self.shorten()

    def extract_galaxy_name(self, file_path, prefix, suffix):
        return file_path.split(prefix)[-1].split(suffix)[0]

    def match_files(self, mag):
        """
        Match all files
        :param mag: quasar mag
        :return:
        """
        # [1:] to remove the F
        bic_pattern = f"{self.direc}/all_out_files_bic_{self.filter_type[1:]}_{mag}_full_quasar_fits_files_6p5_exp_{self.exp_time}_FOV_25.0_samp_{self.ss}/BIC_out_*_db.fits"
        db_pattern = f"{self.direc}/all_out_files_{self.filter_type[1:]}_{mag}_full_quasar_fits_files_6p5_exp_{self.exp_time}_FOV_25.0_samp_{self.ss}/quasar_out_*_db.fits"
        db_woq_pattern = f"{self.direc}/all_out_files_{self.filter_type[1:]}_wo_quasar_fits_files_6p5_exp_{self.exp_time}_FOV_25.0_samp_{self.ss}/without_quasar_out_*_db.fits"
        bic_files = sorted(glob.glob(bic_pattern))
        db_files = sorted(glob.glob(db_pattern))
        db_files_woq = sorted(glob.glob(db_woq_pattern))
        print("Matching files with this pattern...")
        print(db_woq_pattern)
        print("and")
        print(db_pattern)
        # Match files based on galaxy names
        if self.photometric_only:
            matched_files = db_files_woq
        else:
            # for db_file in db_files:
            #     name1 = self.extract_galaxy_name(db_file, 'quasar_out_', '_db.fits')
            #     for db_woq_file in db_files_woq:
            #         name2 = self.extract_galaxy_name(db_woq_file, 'without_quasar_out_', '_db.fits')
            #         print(f"Trying to match: {name1} <-> {name2}")
            matched_files = [
                (db_file, db_woq_file)
                for db_file in db_files
                for db_woq_file in db_files_woq
                if self.extract_galaxy_name(db_file, 'quasar_out_', '_db.fits') == self.extract_galaxy_name(db_woq_file, 'without_quasar_out_', '_db.fits')
            ]
        return matched_files, bic_files

    def process_quasar(self):
        mag = self.quasar_mag_dict[self.quasar_name_long]["quasar"][self.filter_type][0]
        mag = str(mag).replace(".", "p")

        file_to_check = f"{self.psf_mag}psf_collated_data/truth_mag_wq_{quasar_name}_{self.filter_type}_{self.exp_time}.npy"
        print("File to check")
        print(file_to_check)
        if os.path.exists(file_to_check) or self.photometric_only:
            print("Files exist.")
            self.keys_to_load = [
                f"truth_mag_wq_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}",
                f"truth_mag_wq_likelihoods_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}",
                f"truth_mag_sigs_wq_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}",
                f"indices_used_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}"
            ]
            if self.bic:
                self.keys_to_load.append(f"BIC_DIFF_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}")
            if not self.photometric_only:
                self.keys_to_load.extend([f"observed_mag_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}",
                f"observed_mag_likelihoods_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}",
                f"observed_mag_sigs_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}"])
            self.load_results()
        else:
            print("Files do not exist.")

            matched_files, bic_files = self.match_files(mag)
            print("MATCHED FILES")
            print(len(matched_files))
            print("FULL LENGTH OF INDICES")
            print(self.indices_full_length)

            # Define key prefix for storing results
            key_prefix = f"{quasar_name}_{self.filter_type}_{self.exp_time}"

            # Process matched files
            for matched_file in matched_files:
                if self.photometric_only:
                    db_woq_file = matched_file
                    gal_ind_wq = self.extract_galaxy_name(db_woq_file, 'without_quasar_out_', '_db.fits')

                else:
                    db_file, db_woq_file = matched_file
                    gal_ind_wq = self.extract_galaxy_name(db_woq_file, 'without_quasar_out_', '_db.fits')
                    gal_ind = self.extract_galaxy_name(db_file, 'quasar_out_', '_db.fits')
                    assert gal_ind == gal_ind_wq

                try:
                    # Load databases for the galaxy
                    hdul_wo_quasar = fits.open(db_woq_file)
                    db_woq = load_database(db_woq_file)
                    if not self.photometric_only:
                        db = load_database(db_file)
                        hdul_quasar = fits.open(db_file)

                    if self.bic and bic_files:
                        # Load fits files if BIC is required
                        bic_file = bic_files[0] if bic_files else None
                        hdul_bic = fits.open(bic_file)



                    if self.bic and bic_file:
                        ln_bic = hdul_bic[1].data["lnprobability"]
                        ln_quasar = hdul_quasar[1].data["lnprobability"]
                        ln_wo_quasar = hdul_wo_quasar[1].data["lnprobability"]

                        bic_wo_quasar = 1 * np.log(self.num_samples) - 2 * np.min(ln_bic)
                        bic_quasar = 8 * np.log(self.num_samples) - 2 * np.min(ln_quasar)

                        diff = bic_wo_quasar - bic_quasar  # Save or use as needed

                except Exception as e:
                    print(
                        f"-----------------\nError processing galaxy {gal_ind}: {e}\nDB file: {db_file}\n-----------------")
                    continue  # Continue with the next galaxy



                if self.bic:
                    self.results[f"BIC_{key_prefix}"].append(diff)

                if not self.photometric_only:
                    # Extract magnitude traces WTIH QUASAR
                    trace_name_obs = "2_Sersic_mag"
                    trace_obs = _get_trace(trace_name_obs, db)
                    self.results[f"observed_mag_{key_prefix}"].append(np.median(trace_obs))
                    # print(np.median(trace_obs))
                    self.results[f"observed_mag_likelihoods_{key_prefix}"].append(trace_obs.data)
                    self.results[f"observed_mag_sigs_{key_prefix}"].append(np.std(trace_obs))

                    trace_name_true = "1_Sersic_mag"
                    trace_true = _get_trace(trace_name_true, db_woq)
                    self.results[f"truth_mag_wq_{key_prefix}"].append(np.median(trace_true))
                    self.results[f"truth_mag_wq_likelihoods_{key_prefix}"].append(trace_true.data)
                    self.results[f"truth_mag_sigs_wq_{key_prefix}"].append(np.std(trace_true))
                    self.results[f"indices_used_{key_prefix}"].append(gal_ind_wq)

            # Update final indices and save results
            indices_used = self.results[f"indices_used_{key_prefix}"]
            self.indices_used = np.array(indices_used, dtype=int)
            print(f"Found {len(matched_files)} files for {key_prefix}.")
            self.save_results()

    def run_analysis(self):
        self.process_quasar()
        self.shorten() # this doesn't matter anymore as the larger and 108k version are exactly the same indices as long as they're sorted
        try:
            self.create_plots()
        except Exception as e:
            print(e)

    def save_results(self):
        # Save results to .npy files
        for key, value in self.results.items():
            print("Saving ", f"{self.psf_mag}psf_collated_data/{key}.npy")
            np.save(f"{self.psf_mag}psf_collated_data/{key}.npy", value)

    def load_results(self):
        # Load results from .npy files
        self.results = {}
        for key in self.keys_to_load:  # Assuming you have a list of keys to load
            print("Loading ", f"{self.psf_mag}psf_collated_data/{key}.npy")
            self.results[key] = np.load(f"{self.psf_mag}psf_collated_data/{key}.npy", allow_pickle=True)
        self.indices_used = np.asarray(self.results[f"indices_used_{self.quasar_name_long}_{self.filter_type}_{self.exp_time}"], dtype=int)

    def shorten(self):
    # make sure only including indices <108k
        print("self.indices_used")
        print(self.indices_used)
        if np.max(self.indices_used) > 108000:
            mask_less_than_108k = self.indices_used < 108000
            for key, value in self.results.items():
                self.results[key] = np.asarray(self.results[key])[mask_less_than_108k]
            self.indices_used = self.indices_used[mask_less_than_108k]

    def create_plots(self):
        # Create and save histogram plots of observed magnitudes
        plt.figure()
        for i in range(len(self.results[f"truth_mag_wq_{self.quasar_name_long}_{self.filter_type}"])):
            plt.hist(self.results[f"observed_mag_likelihoods_{self.quasar_name_long}_{self.filter_type}"][i], label=f"{self.results[f'truth_mag_wq_{self.quasar_name_long}_{self.filter_type}'][i]:.2f}", density=True)
        plt.legend()
        plt.xlabel("Observed Magnitudes")
        plt.ylabel("Probability")
        plt.savefig(f"../image_pngs/bins_all_quasars.png")
        plt.close()

        # Generate scatter plot of observed vs. true galaxy magnitudes
        fig, ax1 = plt.subplots()
        sc = ax1.scatter(self.results[f"truth_mag_wq_{self.quasar_name_long}_{self.filter_type}"], self.results[f"observed_mag_{self.quasar_name_long}_{self.filter_type}"], c=np.asarray(self.results[f"truth_mag_wq_{self.quasar_name_long}_{self.filter_type}"]) - np.asarray(self.results[f"observed_mag_{self.quasar_name_long}_{self.filter_type}"]))
        ax1.set_xlabel(r'$\mathbf{m_{truth}}$')
        ax1.set_ylabel(r'$\mathbf{m_{\rm obs}}$')
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
        cbar.set_label(r'Residuals ($\mathbf{m_{truth}} - \mathbf{m_{obs}}$)')
        plt.savefig(f"../posteriors/residuals_{self.quasar_name_long}_{self.filter_type}.png", dpi=500)
        plt.close()

        self.make_photometric_plot_appendix()

    def make_photometric_plot_appendix(self, save_residuals=False):
        # Construct key and file paths
        key_prefix = f"{self.quasar_name_long}_{self.filter_type}_{self.exp_time}"
        photo_file = f"/fred/oz183/sberger/paper_1_bluetides/main_scripts/width_by_rad_{int(z)}_JWST.NIRCAM.{self.filter_type}_mag_all_order_luminous_BH.npy"

        # Ensure indices are integers
        self.indices_used = np.array(self.indices_used).astype(int)
        print("length")
        print(len(self.indices_used))

        # Load photometric mags
        photo_mags = np.load(photo_file).flatten()[self.indices_used]
        truth_mags = np.asarray(self.results[f"truth_mag_wq_{key_prefix}"])
        residuals = truth_mags - photo_mags

        # Save large residuals if needed
        if save_residuals:
            indices = np.where(np.abs(residuals) > 1)[0]
            residual_vals = residuals[indices]
            np.save(f"residual_vals_truth-photo_{key_prefix}.npy", residual_vals)

        # Plot
        fig, ax1 = plt.subplots()
        min_mag = np.min(truth_mags)
        max_mag = np.max(truth_mags)
        one_to_one = np.linspace(min_mag, max_mag, 100)

        ax1.plot(one_to_one, one_to_one, color="black", linestyle="--", label="1-1")
        sc = ax1.scatter(truth_mags, photo_mags, c=residuals)
        ax1.set_xlabel(r'$\mathbf{\rm m_{truth}}$', fontsize=18)
        ax1.set_ylabel(r'$\mathbf{\rm m_{photo}}$', fontsize=18)
        ax1.tick_params(axis='both', labelsize=14)
        ax1.invert_xaxis()
        ax1.invert_yaxis()
        xlim = ax1.get_xlim()

        plt.title(f"Galaxies used for {self.quasar_name_long} in {self.filter_type}", fontsize=16)
        cbar = plt.colorbar(sc, ax=ax1, orientation='vertical')
        cbar.set_label(r'Residuals ($\mathbf{m_{\rm truth}} - \mathbf{m_{\rm photo}} $)', fontsize=16)
        cbar.ax.tick_params(labelsize=14)

        plt.legend(fontsize=13)
        plt.xlim(xlim[0], None)
        plt.savefig(f"photometric_residuals_{self.quasar_name_long}_{self.filter_type}.png", dpi=200)
        plt.close()
        print("Made photometric vs Sersic residual plot...")


def get_wavelength(filter_name):
    wavelengths = {
        'F115W': 1.15,
        'F150W': 1.50,
        'F200W': 2.00,
        'F356W': 3.56
    }
    return wavelengths.get(filter_name, 0)



# Example usage
if __name__ == "__main__":
    z = 6.560000154967445
    import quasar_dict

    ### getting photometric plot for where we have the most sampled galaxies (19.5 quasar mag)
    # These directories have the most files
    # 17109 (FIVE FILES PER HOST) all_out_files_356W_wo_quasar_fits_files_6p5_exp_10000_FOV_25.0_samp_2/
    # 10768 (FIVE FILES PER HOST) all_out_files_356W_19p5_full_quasar_fits_files_6p5_exp_10000_FOV_25.0_samp_2/
    # quasar_test_mag = 19.5
    # quasar_data = {
    #     "test": {
    #         "quasar": {
    #             "356W": (quasar_test_mag, (None,))
    #         },
    #         "galaxy": {
    #             "356W": (None, (None,))
    #         },
    #         "category": "ding"
    #
    #     }
    # }
    # #
    # analysis = QuasarParser("test", quasar_data, exp_time=10000, filter_type="356W", ss=2, bic=False, photometric_only=True)
    # # analysis.run_analysis()
    # analysis.make_photometric_plot_appendix(save_residuals=True)
    # exit()
    ### Below is the normal extraction of parameters, above is the huge 19.5 run
    quasar_data = quasar_dict.quasar_data
    exposure_dict = quasar_dict.exposure_times_dict
    for quasar_name, data in quasar_data.items():
        category = data["category"].lower()
        if quasar_name != "J2236+0032":
            continue
        # if "J0844-01" not in quasar_name:
        #     continue
        # Iterate over all filters in the 'quasar' section for this object
        for filter_type in data["quasar"]:
            if filter_type != "F356W":
                continue
            print(f"TRYING QUASAR {quasar_name} in {filter_type}")
            exp_time = exposure_dict[category][filter_type]
            if exp_time is None:
                print(f"Warning: No exposure time found for quasar {quasar_name} using filter {filter_type}")
                continue
            wavelength = get_wavelength(filter_type)
            ss_value = 2 if wavelength > 3 else 1
            print(ss_value)
            analysis = QuasarParser(quasar_name, quasar_data, exp_time=exp_time, filter_type=filter_type, ss=ss_value, photometric_only=False)
            analysis.run_analysis()
            analysis.make_photometric_plot_appendix(save_residuals=False)
