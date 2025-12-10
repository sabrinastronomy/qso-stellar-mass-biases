import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import SynthObs
import SynthObs.Morph
from SynthObs.SED import models
from SynthObs.Morph import measure
from SynthObs.Morph import images
from SynthObs.Morph import PSF
import FLARE.filters
from matplotlib.patches import Circle
import pandas as pd
from synphot import etau_madau
from mpl_toolkits.axes_grid1 import make_axes_locatable
import make_background
from photutils import aperture_photometry
from photutils import CircularAperture
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.modeling.functional_models import Gaussian2D

def create_circular_mask(h, w, center=None, radius=None):
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def ABtonJy(mAB):
    return 10 ** (-0.4 * (mAB - 31.4))

def plot_PSF(LPSF, super_samp, f, noise, background):
    LPSF_mag = LPSF
    LPSF = ABtonJy(LPSF)
    print(f"PSF has mag {LPSF_mag} and flux {LPSF} nJy.")
    filt_str = (f.split('.')[-1])

    img = images.observed(f, cosmo, z, target_width_arcsec=width, smoothing=False, verbose=True, PSF=PSFs[f],
                          super_sampling=super_samp).particle(np.array([0.]), np.array([0.]), np.array([LPSF]),
                                                              centre=[0, 0, 0]) # OUTPUT IN nJy

    if background:
        # create background image object (cutoutwidth in pixels)
        ap_sig = 10 # 10 sigma observation
        background_object = make_background.Background(zeropoint=zeropoint, pixel_scale=pixel_scale / super_samp,
                                                       aperture_f_limit=aperture_f_limit, aperture_significance=ap_sig,
                                                       aperture_radius=aperture_radius, verbose=True)
        img_bkg = background_object.create_background_image(Npixels * super_samp)



        img_bkg_data = img_bkg.bkg * nJy_to_es
        bkg_sigma = background_object.pixel.noise_es * np.ones_like(img_bkg.bkg)

    mask = create_circular_mask(Npixels * super_samp, Npixels * super_samp, center=None,
                                radius=np.floor(Npixels * super_samp / 2))
    img_data = img.super.data * nJy_to_es


    if noise:
        # Add shot noise to full noise
        full_img = img_data * exp_time
        full_img[full_img < 0] = 0
        noisy_full_img = np.random.poisson(full_img)
        img_data = noisy_full_img / exp_time

    if background:
        y, x = np.mgrid[0:len(img_bkg.bkg), 0:len(img_bkg.bkg)]
        gauss = Gaussian2D(np.max(img_data) / 5000, len(img_bkg.bkg) / 2 - 1.5, len(img_bkg.bkg) / 2 - 1.5, 1.5, 1.5)(x,
                                                                                                                      y)
        print('Center loc: ', len(img_bkg.bkg) / 2 - 1.5)
        ivm = 1 / ((bkg_sigma * super_samp) ** 2 + (gauss))
    else: # IVM just 0s
        ivm = np.full_like(img_data, 0.01)


    hdu = fits.PrimaryHDU(mask * img_data) # saving main psf
    hdu_ivm = fits.PrimaryHDU(mask * ivm) # saving IVM for psfmc


    if background and noise:
        hdu.writeto(
            '/fred/oz183/sberger/paper_2_obs_bias/psf_fits_files/sci_PSF_JWST_{}_SN_{}ss_{}s_mag_{}.fits'.format(filt_str,
                                                                                                     super_samp,
                                                                                                    exp_time, LPSF_mag),
            overwrite=True)
        hdu_ivm.writeto(
            '/fred/oz183/sberger/paper_2_obs_bias/psf_fits_files/ivm_PSF_JWST_{}_SN_{}ss_{}s_mag_{}.fits'.format(filt_str,
                                                                                                    super_samp,
                                                                                                    exp_time, LPSF_mag),
            overwrite=True)
        print("Saved PSF...")
        print('/fred/oz183/sberger/paper_2_obs_bias/psf_fits_files/sci_PSF_JWST_{}_SN_{}ss_{}s_mag_{}.fits'.format(filt_str,
                                                                                                     super_samp,
                                                                                                    exp_time, LPSF_mag))
    else: # NOISELESS PSF
        psf_name = '/NOISELESS_sci_PSF_JWST_{}_SN_{}ss_{}s.fits'.format(filt_str, super_samp, exp_time)
        ivm_psf_name = '/fred/oz183/sberger/paper_2_obs_bias/psf_fits_files/NOISELESS_ivm_PSF_JWST_{}_SN_{}ss_{}s.fits'.format(filt_str, super_samp, exp_time)

        print(psf_name)
        print(ivm_psf_name)
        hdu.writeto(psf_name,
            overwrite=True)
        hdu_ivm.writeto(ivm_psf_name,
            overwrite=True)
        print("Noiseless PSF created!")

    return


if __name__ == '__main__':
    # Setup
    cosmo = FLARE.default_cosmo()
    z = 6.560000154967445 # exact BlueTides redshift
    if '--noise' in sys.argv:
        noise = True
    else:
        noise = False

    if '--background' in sys.argv:
        background = True
    else:
        background = False

    filter = sys.argv[1]
    filter_str = str(filter)
    super_samp = int(sys.argv[2])
    exp_time = int(sys.argv[3]) #3100  # 10000

    model = models.define_model('BPASSv2.2.1.binary/ModSalpeter_300') # DEFINE SED GRID, I don't think this matters since we're just grabbing wavelengths

    F = FLARE.filters.add_filters([filter_str], new_lam=model.lam * (1. + z))
    PSFs = PSF.Webb([filter_str], resampling_factor=5)  # creates a dictionary of instances of the webbPSF class

    width = 3  # 8.33 #size of cutout in ''  #MUST RESULT IN EVEN NUMBER OF PIXELS
    FOV = width / cosmo.arcsec_per_kpc_proper(z).value  # size of cutout in kpc
    smoothing = None  # ('adaptive',60)

    pixel_scale = FLARE.filters.pixel_scale[filter_str]  # arcsec/pixel (for NIRCam SW)
    Npixels = int(width / pixel_scale)  # 20#width of image / resolution
    # background setup
    aperture_radius = 2.5 * pixel_scale  # aperture radius in arcsec
    zeropoint = 25.946
    nJy_to_es = 1E-9 * 10 ** (0.4 * (zeropoint - 8.9))

    ## JWST noise estimates are taken from here:
    # https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-sensitivity
    ### THESE ARE ALL IN NJY ####
    aperture_flux_limits_12p3ks = {'JWST.NIRCAM.F115W': 10.496, 'JWST.NIRCAM.F250M': 26.555}
    aperture_flux_limits_10ks = {'JWST.NIRCAM.F150W': 9.267, 'JWST.NIRCAM.F356W': 12.02,
                                      "JWST.NIRCAM.F200W": 8.215}
    # self.aperture_flux_limits_8p2ks = {'JWST.NIRCAM.F115W': 12.57 , 'JWST.NIRCAM.F250M': 32.747} # not used in new proposal
    aperture_flux_limits_5p9ks = {'JWST.NIRCAM.F200W': 10.6937}

    aperture_flux_limits_5p8ks = {'JWST.NIRCAM.F200W': 11.186, 'JWST.NIRCAM.F444W': 21.227}
    aperture_flux_limits_4p3ks = {'JWST.NIRCAM.F115W': 17.6375}  # exactly 4381s

    aperture_flux_limits_4p1ks = {'JWST.NIRCAM.F115W': 18.06, 'JWST.NIRCAM.F200W': 13.005,
                                       'JWST.NIRCAM.F250M': 47.105, 'JWST.NIRCAM.F480M': 69.573}
    aperture_flux_limits_3p1ks = {'JWST.NIRCAM.F150W': 17.15, 'JWST.NIRCAM.F356W': 21.73}  # exactly 3100s
    aperture_flux_limits_1p5ks = {'JWST.NIRCAM.F356W': 31.45551}  # exactly 1578s

    for psf_mag in np.arange(17.5, 22, 0.5):
    # for exp_time in [1578, 3100, 10000]:
        if exp_time == 12300:
            aperture_f_limit = aperture_flux_limits_12p3ks[filter_str]
        elif exp_time == 10000:
            aperture_f_limit = aperture_flux_limits_10ks[filter_str]
        elif exp_time == 8200:
            aperture_f_limit = aperture_flux_limits_8p2ks[filter_str]
        elif exp_time == 5800:
            aperture_f_limit = aperture_flux_limits_5p8ks[filter_str]
        elif exp_time == 5959:
            aperture_f_limit = aperture_flux_limits_5p9ks[filter_str]
        elif exp_time == 4381:
            aperture_f_limit = aperture_flux_limits_4p3ks[filter_str]
        elif exp_time == 4100:
            aperture_f_limit = aperture_flux_limits_4p1ks[filter_str]
        elif exp_time == 3100:
            aperture_f_limit = aperture_flux_limits_3p1ks[filter_str]
        elif exp_time == 1578:
            aperture_f_limit = aperture_flux_limits_1p5ks[filter_str]
        else:
            print("Exposure time not yet supported.")
            exit()
        print('Aperture flux limit ', aperture_f_limit, 'Exp time ', exp_time)
        FPSF = psf_mag  # MAG AB! --> this corresponds to about 19.5 m_ab (6e4 nJy) and 17 m_ab (6e5 nJy), before you used 6e5
        plot_PSF(FPSF, super_samp, filter_str, noise, background)


