from psfMC import model_galaxy_mcmc
import os
import glob
from mpi4py import MPI
from multiprocessing import Pool
import numpy as np
import sys
from config_file_dict_maker import QuasarSimulationConfig
import argparse


# Set up argument parser
parser = argparse.ArgumentParser(description="Run MCMC fitting on a model file.")
parser.add_argument(
    "model_file",
    type=str,
    help="Path to the model file (Python script) to be processed."
)

args = parser.parse_args()
model_file = args.model_file

# These are additional parameters for the MCMC fitter. Number of iterations,
# number of burn-in iterations (which are discarded)
_mcparams = {'burn': 1000, 'iterations': 2000, 'chains': 128}
# _mcparams = {'burn': 500, 'iterations': 1000, 'chains': 128}
# _mcparams = {'burn': 200, 'iterations': 400, 'chains': 128} # quick run


output_name = model_file.replace('model', 'out').replace('.py', '')
model_galaxy_mcmc(model_file, output_name=output_name, **_mcparams)
