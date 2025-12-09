"""
This script is used to check the images made by psfmc and SynthObs in a panel plot.
"""
from astropy.io import fits
from psfMC import model_galaxy_mcmc, load_database
# from psfMC.analysis import plot_hist, corner_plot
from psfMC.analysis.plotting import _get_trace
import glob
import matplotlib

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, SymLogNorm

import subprocess
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import argparse

import FLARE.filters
import FLARE
from SynthObs.SED import models
import importlib
import numpy as np

from astropy.modeling.models import Sersic2D

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 16})

zp = 25.9463

import os
import numpy as np
import matplotlib.pyplot as plt
import corner
from warnings import warn
from matplotlib.ticker import ScalarFormatter, MaxNLocator

import numpy as np
import matplotlib.pyplot as plt
import corner
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import os

import numpy as np
import matplotlib.pyplot as plt
import corner
from matplotlib.ticker import ScalarFormatter, MaxNLocator
import os


def corner_plot_w_quasar(db_file, disp_parameters=None, save=False,
                                  skip_zero_variance=True, filter_walkers=10, **kwargs):
    """
    Make an improved corner plot of all stochastically sampled parameters.
    """

    db = load_database(db_file)
    disp_name = os.path.splitext(os.path.basename(db_file))[0]

    # Improved label list with better formatting
    full_labels = [
        r"Point Source $x$ [pix]",
        r"Point Source $y$ [pix]",
        r"Point Source mag",
        r"Sérsic PA [deg]",
        r"Sérsic mag",
        r"Sérsic ind",
        r"Sérsic $R_{\mathrm{ea}}$ [pix]",
        r"Sérsic $R_{\mathrm{eb}}$ [pix]",
        r"Sérsic $x$ [pix]",
        r"Sérsic $y$ [pix]",
        r"Samples"
    ]

    available_col_names = db.colnames

    if disp_parameters is None:
        display_col_names = list(available_col_names)
        for col in ['lnprobability', 'walker']:
            if col in display_col_names:
                display_col_names.remove(col)
    else:
        missing = set(disp_parameters) - set(available_col_names)
        if missing:
            raise ValueError(f'Missing traces: {missing}')
        display_col_names = list(disp_parameters)

    # Get traces
    traces = []
    for name in display_col_names:
        trace = _get_trace(name, db)
        traces.append(trace)

    flat_traces = np.column_stack(traces)
    n_dim = len(full_labels)

    print(f"flat_traces shape: {flat_traces.shape}")
    print(f"Number of parameters: {n_dim}")

    # Set up the plot with better styling
    plt.style.use('default')  # Start with clean style

    # Create corner plot with improved settings
    fig = corner.corner(
        flat_traces,
        labels=full_labels,  # Use full labels
        max_n_ticks=1,  # More ticks for better readability
        label_kwargs={"fontsize": 14, "labelpad": 15},
        title_kwargs={"fontsize": 12, "pad": 10},
        range=[0.99] * n_dim,
        figsize=(min(3 * n_dim, 16), min(3 * n_dim, 16)),  # Dynamic sizing with max limit
        bins=30,  # More bins for smoother histograms
        smooth=True,  # Smooth the 2D contours
        plot_density=True,  # Show density instead of just contours
        plot_contours=True,
        fill_contours=True,
        color='steelblue',
        hist_kwargs={
            'color': 'steelblue',
            'alpha': 0.7,
            'edgecolor': 'black',
            'linewidth': 0.5
        },
        contour_kwargs={
            'colors': ['steelblue'],
            'linewidths': 1.0
        },
        contourf_kwargs={
            'colors': ['lightsteelblue', 'steelblue'],
            'alpha': 0.6
        },
        **kwargs
    )

    # Get axes array for easier manipulation
    axes = np.array(fig.axes).reshape((n_dim, n_dim))

    # Force no scientific notation on ALL axes - apply after corner plot creation
    for ax in fig.get_axes():
        # Create formatter that absolutely prevents scientific notation
        formatter = ScalarFormatter(useMathText=False)
        formatter.set_scientific(False)
        formatter.set_useOffset(False)
        formatter.set_powerlimits((-99, 99))

        # Apply to both axes
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)

        # Force refresh of tick labels
        ax.ticklabel_format(style='plain', axis='both')

    # Improve formatting for each subplot
    for i in range(n_dim):
        for j in range(n_dim):
            ax = axes[i, j]

            # Skip upper triangle
            if j > i:
                ax.set_visible(False)
                continue

            # Format diagonal (histogram) plots
            if i == j:
                # Clean up histogram appearance
                ax.tick_params(
                    axis='both',
                    which='major',
                    labelsize=9,
                    direction='inout',
                    length=4,
                    width=1,
                    pad=6  # Add padding between ticks and labels
                )
                ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

                # Remove y-axis labels and ticks for histograms
                ax.set_ylabel('')
                ax.tick_params(labelleft=False, left=False)

                # Limit number of ticks for histograms
                ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))

            # Format off-diagonal (scatter) plots
            else:
                ax.tick_params(
                    axis='both',
                    which='major',
                    labelsize=9,
                    direction='inout',
                    length=4,
                    width=1,
                    pad=30  # Add padding between ticks and labels
                )
                ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

                # Limit number of ticks
                ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))

            # Only show x-labels on bottom row
            if i != n_dim - 1:
                ax.set_xlabel('')
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(ax.get_xlabel(), fontsize=11, labelpad=15)

            # Only show y-labels on left column (and not for diagonal)
            if j != 0 or i == j:
                ax.set_ylabel('')
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(ax.get_ylabel(), fontsize=11, labelpad=15)

    # Add a main title
    fig.suptitle('Parameter Posterior Distributions', fontsize=16, y=0.98)

    # Improve overall layout with more spacing
    plt.tight_layout()
    plt.subplots_adjust(top=0.94, bottom=0.12, left=0.12, right=0.98, hspace=0.15, wspace=0.12)

    if save:
        filename = f"{disp_name}_corner_improved.pdf"
        fig.savefig(filename, bbox_inches="tight", dpi=300)
        print(f"Plot saved as {filename}")
    else:
        plt.show()

    return fig


def get_trace_values(trace_name, db, image_data, pkpc_side_fov=25):
    trace = _get_trace(trace_name, db)
    return np.average(trace)

def get_sersic_trace_values(db):
    print(db)
    try:
        amplitude = np.median(_get_trace("2_Sersic_mag", db))
        reff = np.median(_get_trace("2_Sersic_reff", db))
        theta = np.median(_get_trace("2_Sersic_angle", db)) * np.pi / 180 # radians
    except:
        amplitude = np.average(_get_trace("1_Sersic_mag", db))
        reff = np.median(_get_trace("1_Sersic_reff", db))
        theta = np.median(_get_trace("1_Sersic_angle", db)) * np.pi / 180 # radians
    return amplitude, reff, theta

    # if "reff" in trace_name:
    #     side_pixels = np.shape(image_data)[0]
    #     kpc_per_pixel = pkpc_side_fov / side_pixels
    #     print("effective radius [kpc]")
    #     print(np.average(trace) * kpc_per_pixel)
    #     print(np.std(trace) * kpc_per_pixel)
    # if "mag" in trace_name:
    #     print("sersic mag [mag]")
    #     print(np.average(trace))
    #     print(np.std(trace))

def convert_e_s_2_Mjy_sr(flux, curr_filter, both=True): # per pixel
    filter_flare = "JWST.NIRCAM.F" + filter_type.replace("noiseless_", "")
    pixel_scale = FLARE.filters.pixel_scale[filter_flare] # arcsec/pixel (for NIRCam SW)
    print(pixel_scale)
    if curr_filter == "F356W":
        print("ASSUMING SUPERSAMPLING F356W")
        pixel_scale /= 2
    if both:
        flux_MJy = flux / (10**(6 + (0.4*(zp-8.9)))) # e/s to MJy
        pixel_scale_sr = (pixel_scale * (np.pi/(3600*180)))**2
        flux_MJy_sr =  flux_MJy / pixel_scale_sr
        return flux_MJy_sr, pixel_scale
    else:
        return pixel_scale

def plot_four_panel(obs_file_without_quasar, obs_file_with_quasar, direc_fits,
                    direc_fits_wo_quasar, filter_type, quasar_mag, quasar_num, direc):
    print("Making corner plot...")
    db_file_full = direc_fits + "_db.fits"
    corner_plot_w_quasar(db_file_full, save=True)

    db_file_woq = direc_fits_wo_quasar + "_db.fits"

    offset = 25
    new_max = -np.inf
    new_min = np.inf

    # Load image files and convert
    image_files = [
        obs_file_without_quasar,
        obs_file_with_quasar,
        direc_fits + "_convolved_model.fits",
        direc_fits + "_point_source_subtracted.fits",
        direc_fits + "_residual.fits"
    ]
    images = []

    for path in image_files[:-1]:  # All except residual
        try:
            data = fits.getdata(path, ext=0)
        except Exception:
            data = np.full((128, 128), 1)
        data, _ = convert_e_s_2_Mjy_sr(data, filter_type)
        center = data.shape[0] // 2
        data = data[center - offset:center + offset, center - offset:center + offset]
        images.append(data)
        new_max = max(new_max, np.max(data))
        new_min = min(new_min, np.min(data))

    # Residual
    try:
        residual = fits.getdata(image_files[-1], ext=0)
        residual = residual[center - offset:center + offset, center - offset:center + offset]
    except Exception:
        residual = np.full((128, 128), 1)
    residual, _ = convert_e_s_2_Mjy_sr(residual, filter_type)
    images.append(residual)

    # Recovered galaxy - true galaxy residual
    images.append(images[3] - images[0])

    titles = ["true galaxy", "data", "model", "recovered galaxy", "data-model"]
    fig, axes = plt.subplots(1, 5, figsize=(20, 6), sharey=True)

    for i, (img, title) in enumerate(zip(images, titles)):
        sigma = np.std(img)
        tick_font_size = 20
        label_font_size = 20

        if i < 4:
            im = axes[i].imshow(
                img,
                cmap='gist_heat',
                norm=SymLogNorm(linthresh=0.001, vmin=1e-4, vmax=new_max, base=10)
            )
            text_color = "white"
        else:
            im = axes[i].imshow(
                img / sigma,
                cmap='bwr',
                vmin=-6,
                vmax=6
            )
            text_color = "black"

        text_y = 10 if "356" in filter_type else 5
        axes[i].text(1, text_y, title, color=text_color, weight='bold', fontsize=label_font_size)

        # Colorbar below each panel
        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes('bottom', size='5%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=tick_font_size)

        axes[i].set_xticks([])
        if i > 0:
            axes[i].yaxis.set_ticks_position('none')
        axes[i].tick_params(labelsize=tick_font_size)  # Ensures tick font size is large even if some are visible

    plt.tight_layout()
    output_file = f"{quasar_num}_{quasar_mag}_{filter_type}_panel.png"
    plt.savefig(output_file)
    print(f"Saved panel to {output_file}")
    plt.close()

    # Sersic parameter printout
    db_woq = load_database(db_file_woq)
    print("sersic params WO QUASAR:", get_sersic_trace_values(db_woq))

    db_full = load_database(db_file_full)
    print("sersic params W QUASAR:", get_sersic_trace_values(db_full))


def main():
    # plot_title_suffix = f"_{quasar}_{filter_type}_{quasar_num}.png"
    # direc_fits = f'/fred/oz183/sberger/paper_1_bluetides/sab_ding/{quasar}/'
    # direc = f'/fred/oz183/sberger/paper_1_bluetides/mock_im_plots_sixth_round/{quasar}/'
    # direc_initial_images = "/fred/oz183/sberger/paper_1_bluetides/acc_z_fits_files_6p5/"
    # model_file = f'/fred/oz183/sberger/paper_1_bluetides/sab_ding/{quasar}/model_JWST_Ding_{filter_type}.py'
    # obs_file = f'/fred/oz183/sberger/paper_1_bluetides/acc_z_fits_files_6p5/exp_time_{exp_time}_ss_{ss}_not_smoothed_{quasar_num}_{quasar}__with_background_JWST.NIRCAM.{filter_type}.fits'
    # db_file = direc_fits + f"out_JWST_Ding_{filter_type}_db.fits"
    #
    parser = argparse.ArgumentParser(description="Process some directories and files.")

    parser.add_argument('--quasar_num', type=str, required=True, help='quasar_num')
    parser.add_argument('--quasar_mag', type=str, required=True, help='quasar_mag')
    parser.add_argument('--filter_type', type=str, required=True, help='filter used')
    parser.add_argument('--ss', type=str, required=True, help='super sampling')
    parser.add_argument('--exp_time', type=str, required=True, help='exposure time')
    parser.add_argument('--direc', type=str, required=False, help='directory to save images in', default="image_pngs")

    args = parser.parse_args()
    number_filter = args.filter_type.replace("noiseless_", "")

    obs_file_without_quasar = f"/fred/oz183/sberger/paper_2_obs_bias/fit_files_gal_only/{number_filter}_wo_quasar_fits_files_6p5_exp_{args.exp_time}_FOV_25.0_samp_{args.ss}/exp_time_{args.exp_time}_without_quasar_ss_{args.ss}_not_smoothed_{args.quasar_num}__with_background_JWST.NIRCAM.F{number_filter}.fits"
    obs_file_with_quasar = f"/fred/oz183/sberger/paper_2_obs_bias/fits_files_w_quasar/{number_filter}_{args.quasar_mag}_full_quasar_fits_files_6p5_exp_{args.exp_time}_FOV_25.0_samp_{args.ss}/exp_time_{args.exp_time}_ss_{args.ss}_not_smoothed_{args.quasar_num}__with_background_JWST.NIRCAM.F{number_filter}.fits"
    direc_fits = f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_out_files_{args.filter_type}_{args.quasar_mag}_full_quasar_fits_files_6p5_exp_{args.exp_time}_FOV_25.0_samp_{args.ss}/quasar_out_{args.quasar_num}"
    print("DIREC FITS")
    print(direc_fits)
    direc_fits_wo_quasar = f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_out_files_{args.filter_type}_wo_quasar_fits_files_6p5_exp_{args.exp_time}_FOV_25.0_samp_{args.ss}/without_quasar_out_{args.quasar_num}"
    return obs_file_without_quasar, obs_file_with_quasar, direc_fits, direc_fits_wo_quasar, args.filter_type, args.quasar_mag, args.quasar_num, args.direc


if __name__ == '__main__':
    z = 6.5
    obs_file_without_quasar, obs_file_with_quasar, direc_fits, direc_fits_wo_quasar, filter_type, quasar_mag, quasar_num, direc = main()
    ## Observation Setup
    importlib.reload(FLARE.filters)

    cosmo = FLARE.default_cosmo()
    ## ASSUMING IT'S OKAY TO USE THIS SAME MODEL INSTANCE MULTIPLE TIMES
    model = models.define_model('BPASSv2.2.1.binary/ModSalpeter_300')  # DEFINE SED GRID
    filter_flare = "JWST.NIRCAM.F" + filter_type.replace("noiseless_", "")
    all_filters_in_FLARE = FLARE.filters.add_filters([filter_flare], new_lam=model.lam * (1. + z))

    plot_four_panel(obs_file_without_quasar, obs_file_with_quasar, direc_fits, direc_fits_wo_quasar, filter_type, quasar_mag, quasar_num, direc)
