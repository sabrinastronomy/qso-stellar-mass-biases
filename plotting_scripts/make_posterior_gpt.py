import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import norm
import scipy.integrate as spi
from scipy.optimize import curve_fit
from matplotlib import cm
import statsmodels.api as sm
from scipy.interpolate import interp1d
from fastkde import fastKDE
import os
import quasar_dict
import traceback

quasar_data = quasar_dict.quasar_data

def set_plot_defaults(font_size=16):
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    matplotlib.rcParams.update({'font.size': font_size})

class Posterior:
    # Generalized main function
    def __init__(self, filter_type, cdf_use=False, pdf_use=False, quasar_name="", residual_threshold=3, num_bins=None, save_dir=None, num_mags_to_try=2, psf_mag=17,
                 interpolate_fit=False, quasar_gaussian_dict={}, plot_gaussian=False, exposure_time=None):
        self.filter_type = filter_type
        self.quasar_cdf_inferred_parameters = quasar_gaussian_dict

        self.quasar_name_long = quasar_name
        # self.quasar_name = quasar_name.partition("+")[0]

        self.cdf_use = cdf_use
        self.pdf_use = pdf_use
        self.interpolate_fit = interpolate_fit
        self.plot_gaussian = plot_gaussian
        self.exposure_time = exposure_time

        assert self.cdf_use != self.pdf_use # check that both are different boolean values

        if save_dir == None:
            self.save_dir = f"../MAY_2025_posteriors/posteriors_{self.filter_type}_{self.quasar_name_long}"
            if not os.path.exists(self.save_dir):
                os.makedirs(self.save_dir)
                print(f"Directory created at: {self.save_dir}")
            else:
                print(f"Directory already exists at: {self.save_dir}")
        else:
            self.save_dir = save_dir

        # Assign values from the new quasar_data dictionary
        self.mag_quasar = quasar_data[self.quasar_name_long]["quasar"][filter_type][0]  # Magnitude of quasar
        self.m_obs_eiger = quasar_data[self.quasar_name_long]["galaxy"][filter_type][0]  # Magnitude of galaxy
        self.m_obs_eiger_errors = (quasar_data[self.quasar_name_long]["galaxy"][filter_type][1][0],
                              quasar_data[self.quasar_name_long]["galaxy"][filter_type][1][0])  # Error tuple

        self.psf_mag = psf_mag

        # Nicely formatted output
        print(f"\nQuasar Name: {self.quasar_name_long}, Filter: {self.filter_type}, Exp Time: {self.exposure_time}")
        print("=" * 40)
        print(f"{'Variable':<20} {'Value'}")
        print("-" * 40)
        print(f"{'mag_quasar':<20} {self.mag_quasar}")
        print(f"{'m_obs_eiger':<20} {self.m_obs_eiger}")
        print(f"{'m_obs_eiger_errors':<20} {self.m_obs_eiger_errors}")
        print("=" * 40)

        # Load data
        self.photo_mag_arr, self.photo_mag_arr_no_rad, self.obs_mag, self.truth_mag_wq, self.obs_mag_posteriors, self.indices_used = self.load_data()
        print(f"# of samples: {len(self.truth_mag_wq)}")
        print("=" * 40)
        if len(self.truth_mag_wq) < 10:
            return

        residual_check = np.abs(self.obs_mag - self.truth_mag_wq)
        mask_out = residual_check < residual_threshold  # mask out obs mags that are too wrong
        num_true = np.sum(mask_out)
        # print(f"Using {num_true} out of {len(mask_out)}")

        self.obs_mag, self.truth_mag_wq, self.obs_mag_posteriors = self.obs_mag[mask_out], self.truth_mag_wq[mask_out], self.obs_mag_posteriors[mask_out]

        # getting bins from histogram
        plt.close()
        from astropy.stats import bayesian_blocks  # Need Astropy for Bayesian Blocks
        # Compute bin edges using Bayesian Blocks
        # bin_edges = bayesian_blocks(self.truth_mag_wq)
        counts, bin_edges, _ = plt.hist(self.truth_mag_wq, bins=num_bins, edgecolor='black', alpha=0.7)
        print("bin_size")
        print(bin_edges[1]-bin_edges[0])
        plt.xlabel(r"$\rm m_{gal,~truth}$", fontsize=16)
        plt.ylabel("Counts", fontsize=16)
        plt.savefig(f"{self.save_dir}/hist_m_truth_{filter_type}_{quasar_name}.png", dpi=150)
        plt.close()
        self.bins = bin_edges
        self.num_bins = len(bin_edges)
        self.bins_range = (np.round(np.min(self.bins), 2), np.round(np.max(self.bins), 2))
        self.density = np.empty(self.num_bins, dtype=object)  # Array that can hold function objects

        # old way to optimize number of bins
        # if num_bins == None:
        #     n = len(truth_mag_wq)
        #     sigma = np.std(truth_mag_wq)
        #     bin_width = 3.49 * sigma / (n ** (1 / 3))
        #     # Calculate the number of bins
        #     num_bins = int((max_value - min_value) / bin_width)



        # sample_mag = np.linspace(np.min(truth_mag_wq), np.max(truth_mag_wq), num_mags_to_try)
        self.sample_mag = np.linspace(self.m_obs_eiger - 5*self.m_obs_eiger_errors[0], self.m_obs_eiger + 5*self.m_obs_eiger_errors[0], num_mags_to_try) # this goes into integrated posterior

        # Ensure m_obs_eiger is in bins
        if not np.any(np.isclose(self.sample_mag, self.m_obs_eiger, atol=1e-8)):  # atol ensures numerical precision
            self.sample_mag = np.sort(np.append(self.sample_mag, self.m_obs_eiger))  # Add and re-sort

        self.posteriors = np.zeros((self.num_bins - 1))
        self.cdfs = np.zeros((self.num_bins - 1))

        self.residuals_arr = np.zeros_like(self.sample_mag)
        self.mag_truths_bin_centers = (self.bins[1:] + self.bins[:-1]) / 2 # this goes into individual posteriors

        bin_size = self.mag_truths_bin_centers[1] - self.mag_truths_bin_centers[0]

        bin_indices = np.digitize(self.truth_mag_wq, self.bins, right=False) - 1  # subtract 1 to make bins 0-indexed
        self.indices_per_bin = [np.where(bin_indices == i)[0].tolist() for i in range(len(self.bins) - 1)]

        # max mag of posterior to use in residuals plot
        self.max_mags = np.zeros_like(self.sample_mag)


        all_posteriors_raw = np.zeros((len(self.sample_mag), len(self.posteriors)))
        all_posteriors_weighted = np.zeros((len(self.sample_mag), len(self.posteriors)))

        all_cdf_raw = np.zeros((len(self.sample_mag), len(self.posteriors)))
        all_cdf_weighted = np.zeros((len(self.sample_mag), len(self.posteriors)))


        gaussian_mag_arr_pdf = np.zeros_like(self.sample_mag)
        gaussian_mag_arr_cdf = np.zeros_like(self.sample_mag)

        for k, m_obs_eiger_curr in enumerate(self.sample_mag):
            # for k, m_obs_eiger_curr in enumerate([m_obs_eiger]):
            # make posterior
            max_posterior, posterior_truth_mag_given_obs_mag, posterior_not_normalized, cdf_truth_mag_given_obs_mag = self.make_posterior(m_obs_eiger_curr)
            self.max_mags[k] = max_posterior
            # multiplying the posterior values by the value of the gaussian at m_obs
            gaussian_mag_wq_pdf = self.pdf_gaussian(m_obs_eiger_curr, mu=self.m_obs_eiger, sigma=self.m_obs_eiger_errors[0])
            gaussian_mag_wq_cdf = self.cdf_gaussian(m_obs_eiger_curr, mu=self.m_obs_eiger, sigma=self.m_obs_eiger_errors[0])

            gaussian_mag_arr_pdf[k] = gaussian_mag_wq_pdf
            gaussian_mag_arr_cdf[k] = gaussian_mag_wq_cdf

            all_posteriors_raw[k] = posterior_truth_mag_given_obs_mag.copy()
            all_cdf_raw[k] = cdf_truth_mag_given_obs_mag

            # print(f"Multiplying by {gaussian_mag_arr_pdf[k]} and {gaussian_mag_arr_cdf[k]}.")
            posterior_truth_mag_given_obs_mag_weighted = posterior_truth_mag_given_obs_mag * gaussian_mag_wq_pdf
            all_posteriors_weighted[k] = posterior_truth_mag_given_obs_mag_weighted

            cdf_truth_mag_given_obs_mag_weighted =  cdf_truth_mag_given_obs_mag * gaussian_mag_wq_cdf
            all_cdf_weighted[k] = cdf_truth_mag_given_obs_mag_weighted

            self.posteriors += posterior_truth_mag_given_obs_mag_weighted
            self.residuals_arr[k] = max_posterior - m_obs_eiger_curr

            self.cdfs += cdf_truth_mag_given_obs_mag_weighted

        # print("FINAL POSTERIORS WEIGHTED")
        # print(self.posteriors)
        self.posteriors = self.posteriors / np.trapz(self.posteriors, dx=self.mag_truths_bin_centers[1]-self.mag_truths_bin_centers[0]) # renormalizing
        mask_zero = self.posteriors != 0
        self.posteriors = self.posteriors[mask_zero]
        self.mag_truths_bin_centers = self.mag_truths_bin_centers[mask_zero]
        ### DEBUGGING WEIGHTED SUM
        # self.weighted_sum_explanation_plot(all_posteriors_raw, "Raw")
        self.weighted_sum_explanation_plot(all_posteriors_weighted, "Weighted", gaussian_mag_arrs=gaussian_mag_arr_pdf)
        ### weighted sum

        plt.scatter(self.mag_truths_bin_centers, self.posteriors, color="black", label="Corrected with BlueTides")
        plt.xlabel(r"$\rm m_{JWST}$", fontsize=16)
        plt.ylabel(r"$\sum_{i} P(m_{\rm truth} \mid m_{\rm obs}) \cdot N(m_{\rm obs, i}, \sigma_{\rm obs})$", fontsize=16)
        plt.axvline(x=self.m_obs_eiger, color='grey', linestyle='--',
                    label=r'Observed with JWST: $\rm m_{obs, gal, JWST} \pm$' + rf'$ = {self.m_obs_eiger} \pm {self.m_obs_eiger_errors[0]}$')
        plt.axvspan(self.m_obs_eiger - self.m_obs_eiger_errors[0], self.m_obs_eiger + self.m_obs_eiger_errors[1], color='gray', alpha=0.1)
        #### overplot fit
        try:
            if self.interpolate_fit and np.all(self.posteriors != 0):
                interpolation_func = interp1d(self.mag_truths_bin_centers[mask_zero], self.posteriors,
                                              kind='quadratic', fill_value="extrapolate")
                x_smooth = np.linspace(
                    min(min(self.sample_mag), min(self.mag_truths_bin_centers)),
                    max(max(self.sample_mag), max(self.mag_truths_bin_centers)),
                    1000
                )
                y_smooth = interpolation_func(x_smooth)

                plt.plot(x_smooth, y_smooth, color="black")
        except:
            print("Interpolation failed...")

        # plotting fitted and observed posterior
        dx_gauss = 0.01
        posterior_observed = np.arange(np.min(self.mag_truths_bin_centers), np.max(self.mag_truths_bin_centers), dx_gauss)

        try:
            # mu_fit, sigma_fit = self.fit_gaussian(self.mag_truths_bin_centers, self.posteriors)
            # # print("max posterior")
            # # print(np.max(self.posteriors))
            # gaussian_mag_inferred_pdf = self.gaussian_fitted(posterior_observed, mu=mu_fit, sigma=sigma_fit, A=np.max(self.posteriors))
            # if self.plot_gaussian:
            #     plt.plot(posterior_observed, gaussian_mag_inferred_pdf, label=r'Inferred with BlueTides: $\rm m_{obs, gal, BlueTides} \pm$', ls="solid", color="black")

            # gaussian of observation values
            gaussian_mag_true_pdf = self.pdf_gaussian(posterior_observed, mu=self.m_obs_eiger, sigma=self.m_obs_eiger_errors[0])
            # gaussian_normalized = gaussian_mag_true/np.sum(gaussian_mag_true)
            print("gaussian sum")
            print(np.trapz(gaussian_mag_true_pdf, dx=dx_gauss))
            print("posterior sum")
            print(np.trapz(self.posteriors, dx=self.mag_truths_bin_centers[1]-self.mag_truths_bin_centers[0]))
            plt.title(quasar_name + " " + filter_type)
            plt.plot(posterior_observed, gaussian_mag_true_pdf, ls="dashed", color="grey") #label=r"$N(\rm m_{obs, gal, JWST}, \sigma)$"
            plt.legend(loc=1, fontsize=16)
            # plt.xlabel(rf"$\int_{{{min_val}}}^{{{max_val}}} P(m_{{\rm truth}} | m_{{\rm obs}}) \, dm_{{\rm obs}}$")
            plt.savefig(f"{self.save_dir}/weighted_{filter_type}_{quasar_name}.png", dpi=200)
            plt.close()
            print("WEIGHTED POSTERIOR SAVED}")

            ### MAKE CDF
            # Compute the CDF by integrating the PDF
            cdf_values = np.cumsum(self.posteriors)

            # Normalize CDF (ensure last value is exactly 1)
            cdf_values /= cdf_values[-1]
            # Generate x values
            plt.scatter(self.mag_truths_bin_centers, cdf_values, color="black",
                        label="Corrected with BlueTides")
            # plt.plot(posterior_observed, gaussian_mag_true_cdf, label=r"$N(\rm m_{obs, gal, JWST}, \sigma)$",
            #          ls="dashed", color="black")

            if self.interpolate_fit and np.all(self.posteriors != 0):
                from scipy.interpolate import PchipInterpolator
                # interpolation_func = interp1d(self.mag_truths_bin_centers, cdf_values,
                #                               kind='quadratic', fill_value="extrapolate")

                interpolation_func = PchipInterpolator(self.mag_truths_bin_centers, cdf_values)

                x_smooth = np.linspace(
                    min(self.mag_truths_bin_centers),
                    max(self.mag_truths_bin_centers),
                    1000
                )
                y_smooth = interpolation_func(x_smooth)
                # interpolation_func_prob_x = interp1d(cdf_values, self.mag_truths_bin_centers,
                #                               kind='quadratic')

                interpolation_func_prob_x = PchipInterpolator(cdf_values, self.mag_truths_bin_centers)
                # plt.close()
                # plt.plot(cdf_values, interpolation_func_prob_x(cdf_values))
                # plt.title("testing")
                # plt.show()
                mu_cdf = interpolation_func_prob_x(0.5) # get value at middle of distribution
                print(mu_cdf)
                ("------------")
                print("OFFSET")
                print(self.quasar_name_long)
                print(self.filter_type)
                offset = self.m_obs_eiger - mu_cdf
                print(offset)
                # self.quasar_cdf_inferred_parameters[(self.quasar_name_long, self.filter_type)] = (mu_fit, sigma_fit, offset)

                print("------------")
                sigma_top = interpolation_func_prob_x(0.84) - mu_cdf # get top sigma
                sigma_bottom = mu_cdf - interpolation_func_prob_x(0.16)  # bottom sigma
                plt.plot(x_smooth, y_smooth, color="black")

                mu_actual = self.m_obs_eiger
                mu_actual_err = self.m_obs_eiger_errors

                self.quasar_cdf_inferred_parameters[(self.quasar_name_long, self.filter_type)] = (mu_cdf, (sigma_bottom, sigma_top), offset, mu_actual, (mu_actual_err[0], mu_actual_err[1]), bin_size)

            # except:
            #     print("Interpolation failed...")
            if (self.quasar_name_long, self.filter_type) in {
                ("J1120+0641", "F356W"),
                ("J0844-0132", "F150W"),
                ("J0911+0152", "F150W"),
                ("J1146-0005", "F150W"),
            }:
                plt.axhline(y=0.14, color='red', linestyle='--',
                        label=r'$\rm m_{obs, lower~limit, 16th~percentile}$')
                loc_legend = 2
            else:
                plt.axhline(y=0.5, color='black', linestyle='--',
                        label=r'$\rm m_{obs, inferred}$')
                loc_legend = 4
            # plt.axvspan(self.m_obs_eiger - self.m_obs_eiger_errors[0], self.m_obs_eiger + self.m_obs_eiger_errors[1], color='gray', alpha=0.1)
            # plt.plot(posterior_observed, gauss_cdf(posterior_observed), label=r"$N(\rm m_{obs, gal, JWST}, \sigma)$", ls="dashed", color="black")
            # plt.plot(self.sample_mag, gaussian_mag_arr, label=r"$N(\rm m_{obs, gal, JWST}, \sigma)$", ls="dashed", color="black")
            plt.xlabel(r"$\rm m_{JWST}$", fontsize=14)
            plt.ylabel(r"$\rm P(\rm m_{truth} > m_{obs, x} | m_{obs})$", fontsize=14)
            plt.axvline(x=self.m_obs_eiger, color='grey', linestyle='--',
                        label=r'Observed with JWST: $\rm m_{obs, gal, JWST}$' + rf'$ = {self.m_obs_eiger} \pm {self.m_obs_eiger_errors[0]}$')
            plt.axvspan(self.m_obs_eiger - self.m_obs_eiger_errors[0], self.m_obs_eiger + self.m_obs_eiger_errors[1],
                        color='gray', alpha=0.1)

            plt.legend(loc=loc_legend, fontsize=10)
            plt.ylim(-0.1, 1.1)
            plt.savefig(f"{self.save_dir}/weighted_CDF_{filter_type}_{quasar_name}.png", dpi=200)
            plt.close()
            print("WEIGHTED CDF SAVED")

        except Exception as e:
            print("Exception occurred:")
            traceback.print_exc()  # Prints the full traceback
            self.quasar_cdf_inferred_parameters[(self.quasar_name_long, self.filter_type)] = (None, None, None, None, None, None)

    def generate_bin_size_latex_table(self, output_filename="bin_size.txt"):
        """
        Generates a LaTeX table of bin sizes for each (quasar, filter) combo and writes it to a text file.

        Parameters:
            quasar_cdf_inferred_parameters (dict): Dictionary with keys (quasar_name, filter_type) and
                                                   values as tuples ending in bin_size.
            output_filename (str): Output file name for the LaTeX table.
        """
        latex_table = [
            r"\begin{tabular}{l c}",
            r"\hline",
            r"\textbf{Quasar-Filter} & \textbf{Bin Width (mag)} \\",
            r"\hline"
        ]

        for (quasar_name, filter_type), params in self.quasar_cdf_inferred_parameters.items():
            bin_size = params[5]
            if bin_size == None:
                continue
            key = f"{quasar_name} {filter_type}"
            latex_table.append(f"{key} & {bin_size:.2f} \\\\")

        latex_table.append(r"\hline")
        latex_table.append(r"\end{tabular}")

        # Write to file
        with open(output_filename, "w") as f:
            f.write("\n".join(latex_table))

        print(f"LaTeX table saved to {output_filename}")

    def plot_stellar_mass(self, indices_used):
        stellarMass = np.load('/fred/oz183/mmarshal//BlueTides/stellarMass_208.npy')[:108001]
        indices_used = np.asarray(indices_used, dtype=int)
        plt.hist(stellarMass[indices_used])
        plt.ylabel("stellar mass [solar masses]")
        plt.title(self.quasar_name_long)

        plt.savefig(f"{self.save_dir}/stellar_mass_{self.quasar_name_long}.png", dpi=150)
        plt.close()
    # Define the Gaussian function
    def gaussian_fitted(self, x, mu=25.72, sigma=0.5, A=1):
        return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

    def pdf_gaussian(self, x, mu=25.72, sigma=0.5, A=1):
        pdf_val_normal = norm.pdf(x, loc=mu, scale=sigma)  # making numpy do the normalization :)
        return pdf_val_normal

    def cdf_gaussian(self, x, mu=25.72, sigma=0.5):
        cdf_val_normal = norm.cdf(x, loc=mu, scale=sigma)  # making numpy do the normalization :)
        return cdf_val_normal
        # return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    def renormalize_cdf(self, cdf):
        return (cdf - np.min(cdf)) / (np.max(cdf) - np.min(cdf))


    # Function to fit a Gaussian to (x, y) data
    def fit_gaussian(self, x_data, y_data):
        """
        Fit a Gaussian function to given (x, y) data.

        Parameters:
            x_data (array-like): x values.
            y_data (array-like): y values.

        Returns:
            tuple: Estimated parameters (A, mu, sigma).
        """

        # Initial guess: A = max(y), mu = mean of x, sigma = std of x
        initial_guess = [np.mean(x_data), np.std(x_data)]

        # Fit the Gaussian curve
        popt, _ = curve_fit(self.gaussian_fitted, x_data, y_data, p0=initial_guess)
        # print(popt)
        return popt  # Returns (mu, sigma)

    # Function to load data, generalized for different paths and filenames
    def load_data(self):
        photo_mag_arr = ""
        photo_mag_arr_no_rad = ""

        # self.filter_type = self.filter_type.replace('F', '')

        # Save results to .npy files
        obs_mag = np.load(f"{self.psf_mag}psf_collated_data/observed_mag_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.npy") # should already be medial from all_quasar_likelihood_arr.py
        obs_mag_posteriors = np.load(f"{self.psf_mag}psf_collated_data/observed_mag_likelihoods_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.npy", allow_pickle=True)
        # obs_mag_posteriors = np.squeeze(obs_mag_posteriors, axis=-1)
        # obs_mag = np.median(obs_mag_posteriors, axis=1)

        obs_mag_sigs = np.load(f"{self.psf_mag}psf_collated_data/observed_mag_sigs_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.npy")

        # truth_mag_wq_posteriors = np.load(f"{self.psf_mag}psf_collated_data/truth_mag_wq_likelihoods_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.npy", allow_pickle=True)
        # truth_mag_wq_posteriors = np.squeeze(truth_mag_wq_posteriors, axis=-1)
        truth_mag_wq = np.load(f"{self.psf_mag}psf_collated_data/truth_mag_wq_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.npy")  # should already be medial from all_quasar_likelihood_arr.py
        # truth_mag_wq = np.median(truth_mag_wq_posteriors, axis=1)
        truth_mag_wq_sigs = np.load(f"{self.psf_mag}psf_collated_data/truth_mag_sigs_wq_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.npy")
        indices_used = np.load(f"{self.psf_mag}psf_collated_data/indices_used_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.npy")

        return photo_mag_arr, photo_mag_arr_no_rad, obs_mag, truth_mag_wq, obs_mag_posteriors, indices_used

    def weighted_sum_explanation_plot(self, all_posteriors, extra_title, gaussian_mag_arrs=[]):
        # Close all previous plots
        plt.close("all")
        if "weighted" in extra_title.lower():
            max_val = np.max(gaussian_mag_arrs)
        else:
            max_val = 1.7

        # Set up a large panel with multiple subplots
        num_plots = len(self.sample_mag)
        cols = 3  # Number of columns in the panel
        rows = int(np.ceil(num_plots / cols))  # Compute required rows

        fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows), sharex=False, sharey=True)
        axes = axes.flatten()  # Flatten in case of a grid

        # print("length of all posteriors:", len(all_posteriors))

        for j in range(num_plots):
            # if np.all(all_posteriors[j] == 0):
                # continue  # Skip empty posteriors

            ax = axes[j]  # Get subplot for current index

            # Generate colors per point, rather than per plot
            point_colors = cm.Spectral(np.linspace(0, 1, len(self.mag_truths_bin_centers)))

            for i in range(len(self.mag_truths_bin_centers)):
                posterior_value = all_posteriors[j][i]


                # Determine marker based on indices_per_bin length
                length = len(self.indices_per_bin[i])
                marker_type = 'x' if length < 2 else 'o'
                if marker_type == "x":
                    ax.scatter(self.mag_truths_bin_centers[i], posterior_value,
                               marker=marker_type, label="1 galaxy in bin", color=point_colors[i])
                else:
                    ax.scatter(self.mag_truths_bin_centers[i], posterior_value, color=point_colors[i],
                               marker=marker_type)

            # Vertical bin lines and other reference lines
            if "weighted" in extra_title.lower():
                ax.set_title(rf"Sample magnitude #{j + 1} weighted by {gaussian_mag_arrs[j]:.3f}", fontsize=16)
            else:
                ax.set_title(f"Sample magnitude #{j + 1} weighted by {gaussian_mag_arrs[j]:.3f}", fontsize=16)
            ax.vlines(self.bins, ymin=0, ymax=max_val, label='Bins', alpha=0.2, color="black")
            ax.axvline(x=self.sample_mag[j], ymin=0, ymax=max_val, color='black', linestyle='--', label=r"$\rm m_{gal, obs, JWST}=$" + rf" {self.sample_mag[j]:.2f}")
            ax.axvline(x=np.min(self.truth_mag_wq), ymin=0, ymax=max_val, linestyle='dotted', label=r'$\rm m_{gal,~min/max,~truth,~sample}$',
                       color='black')
            ax.axvline(x=np.max(self.truth_mag_wq), ymin=0, ymax=max_val, linestyle='dotted', color='black')
            ax.legend(fontsize=16, loc=1)
            ax.set_xlabel(r"$\rm m_{truth}$", fontsize=16)

        # Compute element-wise sum of all posteriors
        summed_posteriors = np.sum(all_posteriors, axis=0)  # Summing across all samples
        # mask_zero = summed_posteriors != 0 # need to match mag_truths_bin_centers
        # summed_posteriors = summed_posteriors[mask_zero]
        ax = axes[-1]
        # last panel BELOW!
        for i in range(len(self.mag_truths_bin_centers)):
            posterior_value = summed_posteriors[i]

            # if posterior_value == 0:
            #     continue

            marker_type = 's' if len(self.indices_per_bin[i]) < 1 else 'o'

            ax.scatter(self.mag_truths_bin_centers[i], posterior_value, color=point_colors[i],
                       marker=marker_type)
        dx_gauss = 0.01
        posterior_observed_mags = np.arange(np.min(self.mag_truths_bin_centers), np.max(self.mag_truths_bin_centers), dx_gauss)
        gaussian_mag_at_sample_mags = self.pdf_gaussian(posterior_observed_mags, mu=self.m_obs_eiger, sigma=self.m_obs_eiger_errors[0])

        ax.plot(posterior_observed_mags, gaussian_mag_at_sample_mags, label=r"$N(\rm m_{gal, obs, JWST}, \sigma)$", ls="dashed",
                     color="black", zorder=-10, alpha=0.5)
        ax.legend(loc=1)
        ax.vlines(self.bins, ymin=0, ymax=max_val, label='Bins', alpha=0.2, color="black")
        ax.set_title(rf"Element-wise {extra_title} Sum of Posteriors", fontsize=16)
        # Add a global title
        # fig.suptitle(rf"Weighted Sum - {extra_title} Posteriors", fontsize=14)
        ax.set_xlabel(r"$\rm m_{truth}$")

        # Adjust layout
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        # Save figure
        fig.savefig(f"{self.save_dir}/{self.num_bins}_posterior_summed_pre_panel_{extra_title}.pdf")
        plt.close("all")
        plt.close(fig)

    # Function to compute the posterior and plot results
    def make_posterior(self, m_obs):
        """
        Generates normalized (sums to 1) P(m_obs|m_truth) for a given m_obs
        :param m_obs: JWST observed magnitude
        :return:
        integral_value: magnitude at peak of posterior
        posterior_truth_mag_given_obs_mag: normalized P(m_obs|m_truth)
        cdf_val: cdf version
        """
        # if np.min(self.obs_mag_posteriors) > m_obs:
        #     return np.nan, np.nan

        fig_hist, ax_hist = plt.subplots(figsize=(12, 8))  # Wider figure to fit legend
        colors = cm.Spectral(np.linspace(0, 1, len(self.indices_per_bin)))  # get colors for each bin
        posterior = np.zeros_like(self.mag_truths_bin_centers)
        cdf_vals = np.zeros_like(self.mag_truths_bin_centers)

        for j, indices_in_bin in enumerate(self.indices_per_bin):
            # print("length", len(indices_in_bin))
            if len(indices_in_bin) < 1:
                posterior[j] = 0
                cdf_vals[j] = 0
                # print("Bin does not have any indices.")
                continue

            mag_obs_likelihood_posteriors_selected = np.concatenate(
                [arr.flatten() for arr in self.obs_mag_posteriors[indices_in_bin]]
            )
            # mag_obs_likelihood_posteriors_selected = self.obs_mag_posteriors[indices_in_bin].flatten()
            # print(mag_obs_likelihood_posteriors_selected)
            # print(np.shape(mag_obs_likelihood_posteriors_selected))

            min_likelihood, max_likelihood = min(mag_obs_likelihood_posteriors_selected), max(mag_obs_likelihood_posteriors_selected)
            likelihood_spacing_arr = np.linspace(min_likelihood, max_likelihood, num=1000)



            # Fit a normal distribution and compute the likelihood
            if not callable(self.density[j]) and self.pdf_use: # only do if self.density function at in dex(cdf/pdf interpolated function) has NOT been created)
                # Fit the KDE model
                # Compute the KDE using fastKDE
                pdf = fastKDE.pdf(mag_obs_likelihood_posteriors_selected, var_names = ['x'])
                self.density = interp1d(pdf.variable, pdf.values, kind="quadratic", fill_value=0, bounds_error=False)
            elif not callable(self.density[j]) and self.cdf_use:
                self.density[j] = sm.distributions.ECDF(mag_obs_likelihood_posteriors_selected)
            ######################################################
            # do actual posterior calculation
            if self.pdf_use:
                p = self.density[j](m_obs)
                posterior[j] = p

            elif self.cdf_use:
                cdf_val = self.density[j](m_obs)
                if min_likelihood > m_obs or max_likelihood < m_obs:
                    posterior[j] = 0
                    cdf_vals[j] = 0
                    # print("Selected likelihood samples are out of range.")
                    # continue
                if np.isnan(cdf_val):
                    posterior[j] = 0
                    cdf_vals[j] = 0
                    # continue

                dx = 0.005
                derivative_at_p = (self.density[j](m_obs + dx / 2) - self.density[j](m_obs - dx / 2)) / (dx)

                posterior[j] = derivative_at_p
                cdf_vals[j] = cdf_val

            # plot each CDF/PDF
            ax_hist.plot(likelihood_spacing_arr, self.density[j](likelihood_spacing_arr),
                         color=colors[j], label=r"$\rm m_{gal, truth}~=$" +f"{self.mag_truths_bin_centers[j]:.2f} " +  r"($\rm n_{gal}$" + f"= {len(indices_in_bin)})")
        # checking for bad posterior
        # enters this if statement if none of the PDFs/CDFs covered m_obs ever
        # print("posterior_truth_mag_given_obs_mag")
        if np.all(posterior == 0):
            print("Posterior all 0s")
            return np.nan, posterior, posterior, cdf_vals


        ax_hist.axvline(x=m_obs, color='black', linestyle='--', label=f'Observed JWST Galaxy m = {m_obs:.2f}')
        # With these lines:
        ax_hist.legend(fontsize=16, bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()
        ax_hist.set_xlabel(r"All $\rm m_{obs}$ samples in bin", fontsize=16)
        ax_hist.set_ylabel(r"$\rm P(\rm m_{obs,~likelihood~samples} > m_{obs, x} | m_{truth})$", fontsize=16)

        ax_hist.set_title(self.quasar_name_long)

        if self.pdf_use:
            fig_hist.savefig(
                f"{self.save_dir}/pdf_bins4posterior_{m_obs}_{self.bins_range}_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.png",
                dpi=200, bbox_inches='tight')
            plt.close(fig_hist)
        elif self.cdf_use:
            fig_hist.savefig(f"{self.save_dir}/cdf_bins4posterior_{m_obs}_{self.bins_range}_{self.quasar_name_long}_{self.filter_type}_{self.exposure_time}.png",
                dpi=200, bbox_inches='tight')
            plt.close(fig_hist)
        # Plot each point with its corresponding color
        plt.close()
        posterior_truth_mag_given_obs_mag = posterior / np.nansum(posterior)
        posterior_not_normalized = np.copy(posterior)
        for i in range(len(colors)):
            # print(i)
            # print(posterior_truth_mag_given_obs_mag[i])
            length = len(self.indices_per_bin[i])
            if length < 2:
                plt.scatter(self.mag_truths_bin_centers[i], posterior_truth_mag_given_obs_mag[i], color=colors[i], marker='s')

            else:
                plt.scatter(self.mag_truths_bin_centers[i], posterior_truth_mag_given_obs_mag[i], color=colors[i])



        plt.vlines(self.bins, ymin=0, ymax=1, label='Bins', alpha=0.2, color="black")
        plt.axvline(x=m_obs, color='black', linestyle='--', label=f'Observed EIGER Galaxy m = {m_obs}')
        plt.axvline(x=np.min(self.truth_mag_wq), linestyle='dotted', label=r'$\rm m_{gal,~min/max,~truth,~sample}$', color='black')
        plt.axvline(x=np.max(self.truth_mag_wq), linestyle='dotted', color='black')

        if self.interpolate_fit:
            # Interpolation for a smoother curve
            interpolation_func = interp1d(self.mag_truths_bin_centers, posterior_truth_mag_given_obs_mag,
                                          kind='cubic')
            x_smooth = np.linspace(min(self.mag_truths_bin_centers), max(self.mag_truths_bin_centers), 100)
            y_smooth = interpolation_func(x_smooth)
            plt.plot(x_smooth, y_smooth, color="black", alpha=0.9, zorder=-100)

        plt.ylim((None,1))
        plt.title(self.quasar_name_long)
        plt.xlabel(r"$\rm m_{truth}$")
        plt.ylabel(r"$\rm P(m_{truth} | m_{obs, JWST})$")
        plt.legend(loc=0, fontsize=16)
        plt.tight_layout()  # leave space on the right for the legend
        if self.cdf_use:
            print("Saving CDF plot")
            print(f"{self.save_dir}/cdf_posterior_{m_obs}_{self.bins_range}_{self.quasar_name_long}_{filter_type}.png")
            plt.savefig(f"{self.save_dir}/cdf_posterior_{m_obs}_{self.bins_range}_{self.quasar_name_long}_{filter_type}.png", dpi=150)
        elif self.pdf_use:
            plt.savefig(f"{self.save_dir}/pdf_posterior_{m_obs}_{self.bins_range}_{self.quasar_name_long}_{filter_type}.png", dpi=150)
        plt.close()


        ind = np.nanargmax(posterior_truth_mag_given_obs_mag)
        integral_value = self.mag_truths_bin_centers[ind]

        # assert np.isclose(np.sum(posterior_truth_mag_given_obs_mag), 1)
        return integral_value, posterior_truth_mag_given_obs_mag, posterior_not_normalized, cdf_vals

    def generate_latex_table(self):

        table_lines = [
            r"\begin{tabular}{l l c c c}",  # Use raw strings (r"") to avoid escape issues
            r"    \hline",
            r"    \textbf{Quasar} & \textbf{Filter} & \textbf{$\rm m_{obs, gal, JWST}$} & \textbf{$\rm m_{obs, gal, BlueTides}$} & $\rm m_{obs, gal, JWST} - m_{obs, gal, BlueTides}$ \\",
            r"    \hline",
        ]

        for (quasar, filter_band), (mag_inferred, mag_err_inferred, offset, mag_obs, mag_err_obs, _) in self.quasar_cdf_inferred_parameters.items():
            # Formatting magnitude values
            if mag_inferred is not None:
                obs_bluetides = f"${{{mag_inferred:.2f}}}^{{+{mag_err_inferred[1]:.2f}}}_{{-{mag_err_inferred[0]:.2f}}}$"
                obs_jwst = f"${{{mag_obs:.2f}}}^{{+{mag_err_obs[1]:.2f}}}_{{-{mag_err_obs[0]:.2f}}}$"
                offset = "{:.2f}".format(offset)
            else:
                obs_jwst = obs_bluetides = "--"
                offset = "--"

            # Append row to table
            table_lines.append(rf"{quasar} & {filter_band} & {obs_jwst} & {obs_bluetides} & {offset} \\")

        table_lines.append("\hline")
        table_lines.append("\end{tabular}")

        with open("quasar_inferred_params_cdf_latex_table.txt", "w") as file:
            file.write("\n".join(table_lines))
        return table_lines



# Usage example:
if __name__ == "__main__":
    # Set default plot parameters
    set_plot_defaults(font_size=12)
    z = 6.560000154967445
    import quasar_dict
    quasar_data = quasar_dict.quasar_data
    exposure_dict = quasar_dict.exposure_times_dict

    mag_arr_path = ""
    mag_arr_no_rad_path = ""
    quasar_dict = {}
    for quasar_name, data in quasar_data.items():
        print(quasar_name)
        data = quasar_data[quasar_name]
        for filter_type in data["quasar"]:
            exposure_time = exposure_dict[data["category"]][filter_type]
            if (quasar_name, filter_type) not in {
                # ("J1120+0641", "F200W"),
                # ("J159-02"),
                # ("J2255+0251")
                ("J2236+0032", "F356W"),

            }:
                continue
            post = Posterior(filter_type=filter_type, cdf_use=True, pdf_use=False, quasar_name=quasar_name, num_mags_to_try=10, num_bins=20, interpolate_fit=True, quasar_gaussian_dict=quasar_dict, plot_gaussian=False, exposure_time=exposure_time)
            quasar_dict = post.quasar_cdf_inferred_parameters
            post.generate_latex_table()
            post.generate_bin_size_latex_table()