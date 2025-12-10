"""
This script wraps the entire imaging and fitting pipeline.
1) config_file_dict_maker.py
2) get_all_fits_wrapper
3) model_maker.py
4) fit_gal_run.py OR fit_without_gal_run.py
"""
import json
import argparse
import os
import subprocess
from config_file_dict_maker import QuasarSimulationConfig
import numpy as np
import sys

def str_to_bool(s):
    # Function to convert string to boolean written by ChatGPT
    s = s.lower()
    if s in ['true', 't', 'yes', 'y', '1']:
        return True
    elif s in ['false', 'f', 'no', 'n', '0']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Argument parser for mock imaging config file")
    # Required arguments (no defaults provided)
    parser.add_argument('--run_through_psfmc', type=str_to_bool, help='Whether to run entire pipeline (create dictionary, image all galaxies with/without quasar, make model files, fit the psfmc)')
    parser.add_argument('--include_quasar', type=str_to_bool, help='Whether to include quasar')
    parser.add_argument('--filters', help='JWST filters')
    parser.add_argument('--mags_AB', help='Magnitude of quasar(s) in AB system')
    parser.add_argument('--samps', help='Supersampling (usually 1 or 2)')
    parser.add_argument('--exp_time', help='Exposure time [s]')
    parser.add_argument('--indices_to_fit', type=str, help="file name of galaxies (ordered by BHAR) to fit")

    # Optional arguments with defaults
    parser.add_argument('--psf_mag', default=19.5, help='magnitude of PSF to use')
    parser.add_argument('--FOVs', default=25, help='FOV in pkpc (can also be specified as a width instead)')
    parser.add_argument('--redshift', type=str, default=6.5, help='Redshift of BlueTides simulation to use')
    parser.add_argument('--create_background', type=str_to_bool, default=True, help='Create background (default: True)')
    parser.add_argument('--choose_samp_4_me', type=int, default=False, help='Choose sample for me (default: False)')
    parser.add_argument('--smooth', type=str_to_bool, default=False, help='Apply smoothing (default: False)')
    parser.add_argument('--dust', type=str_to_bool, default=True, help='Apply dust effects (default: True)')
    parser.add_argument('--widths', help='Widths parameters (default: None)', default=None)
    parser.add_argument('--use_MPI', type=str_to_bool, default=True, help='Use MPI (default: True)')
    parser.add_argument('--verbose', type=str_to_bool, default=False, help='Verbose output (default: False)')
    parser.add_argument('--rank', type=int, default=None, help='MPI rank (default: None)')
    parser.add_argument('--size', type=int, default=None, help='MPI size (default: None)')
    parser.add_argument('--comm', help='MPI communicator (default: None)', default=None)
    parser.add_argument('--MPI_tasks', type=str, default=128, help="how many MPI tasks")

    parser.add_argument('--shot', type=str_to_bool, default=True, help='Include shot noise (default: True)')
    parser.add_argument('--ding', type=str_to_bool, default=False, help='Ding quasars')
    parser.add_argument('--just_quasar', type=str_to_bool, default=False, help='ONLY IMAGE QUASAR (CAREFUL) no galaxy imaged if this is True')
    parser.add_argument('--noiseless_psf', type=str_to_bool, default=True, help="Do not include noise in your PSF")
    parser.add_argument('--generate_images', type=str_to_bool, default=True, help="Whether or not to image galaxies/quasar and galaxy")
    parser.add_argument('--BIC', type=str_to_bool, default=False, help="BIC, point source ONLY run")
    parser.add_argument('--full', type=str_to_bool, default=False, help="if full=True, wrapper will grab FULL BlueTides rather than 108k subset")

    args = parser.parse_args()
    print("Inputted args parsed:")
    print(args)

    ### putting args in list form that are needed for multiple run capability in make_images_detected_quasar.py
    # i.e., filter, quasar_mag, samp, FOV, width, we take them out of lists in model_maker.py
    # filters
    args.filters = [args.filters]
    # quasar mags
    if args.include_quasar or args.just_quasar:
        args.mags_AB = [float(args.mags_AB)]
    else:
        args.mags_AB = [None]
    # samples
    args.samps = [int(args.samps)]
    # FOVs
    args.FOVs = [float(args.FOVs)]
    # widths
    if args.widths is not None:
        args.widths = [float(args.widths)]
    else:
        args.widths = [None]

    args.exp_time = int(args.exp_time)
    args.psf_mag = str(args.psf_mag)

    args.indices_to_fit = args.indices_to_fit

    # Create a configuration object from the parsed arguments
    return QuasarSimulationConfig(**vars(args))

if __name__ == "__main__":
    config = parse_arguments()

    # Define the output file for the configuration object
    output_file = f"/fred/oz183/sberger/paper_2_obs_bias/src/config_files/config_file_{config.distinguish_string}.json" # saving to config files directory
    # Save the configuration object to a JSON file
    config.save_to_json(output_file)

    print(f"Configuration saved to {output_file}")

    print("Saving psfmc model files and results with this:")
    print(config.distinguish_string)
    #
    if config.generate_images:
        # image all galaxies
        try:
            # Run the subprocess
            fits_process = subprocess.run(
                ["python", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/get_all_fits_wrapper.py",
                 "--json_config_file", output_file], check=True)
            # If check=True, it will raise a CalledProcessError on failure
            print("------------------------------------------------------------------------------------------")
            print("Finished generating images.")
            print("------------------------------------------------------------------------------------------")

        except subprocess.CalledProcessError as e:
            print("------------------------------------------------------------------------------------------")
            print("ERROR generating images.")
            print(e)
            print("------------------------------------------------------------------------------------------")
            sys.exit(1)

    try:
        # generate all model files
        generate_fits = subprocess.run(
            ["python", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/model_maker.py", "--json_config_file",
             output_file], check = True)
        # If check=True, it will raise a CalledProcessError on failure
        print("------------------------------------------------------------------------------------------")
        print("Finished generating model files.")
        print("------------------------------------------------------------------------------------------")
    except subprocess.CalledProcessError as e:
        print("------------------------------------------------------------------------------------------")
        print("ERROR generating model files.")
        print("------------------------------------------------------------------------------------------")
        sys.exit(1)


    if config.run_through_psfmc:
        print("STARTING FITTING")
        # do fitting
        if config.include_quasar:
            fit_galaxy = subprocess.run(["mpirun", "-n", config.MPI_tasks, "python", "-m", "mpi4py", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/fit_gal_run.py", "--json_config_file", output_file], check=True)
        else:
            fit_galaxy = subprocess.run(["mpirun", "-n", config.MPI_tasks, "python", "-m", "mpi4py", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/fit_without_gal_run.py", "--json_config_file", output_file], check=True)
