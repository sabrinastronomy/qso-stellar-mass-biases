from psfMC import model_galaxy_mcmc
import os
import glob
from mpi4py import MPI
from multiprocessing import Pool
import numpy as np
import sys
from config_file_dict_maker import QuasarSimulationConfig
import argparse

"""Script to use psfmc on many galaxies using MPI"""

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(" Running %d parallel MPI processes" % size)

comm.Barrier()
t_start = MPI.Wtime()

# print("group per node")
# print(group_per_node)

parser = argparse.ArgumentParser(description="Which json config file should I use?")
parser.add_argument('--json_config_file', type=str, help='json config file to use')
args = parser.parse_args()
config = QuasarSimulationConfig.load_from_json(args.json_config_file)
filter_type = config.filters
include_quasar = config.include_quasar

filter_first_str = filter_type[0][-4:]
FOVs = config.FOVs[0]  # pkpc
widths = config.widths[0]  # arcseconds
exp_time = config.exp_time
ss = config.samps[0]
exp_time = config.exp_time
indices_to_fit = np.load(config.indices_to_fit)

# These are additional parameters for the MCMC fitter. Number of iterations,
# number of burn-in iterations (which are discarded)
_mcparams = {'burn': 1000, 'iterations': 2000, 'chains': 128}
# _mcparams = {'burn': 500, 'iterations': 1000, 'chains': 128}
# _mcparams = {'burn': 200, 'iterations': 400, 'chains': 128} # quick run

def run_mcmc(model_file):
    output_name = model_file.replace('model', 'out').replace('.py', '')
    if len(glob.glob(output_name +'_*fits*')) > 0:
        print(glob.glob(output_name +'*fits*'))
        print('{} already processed, skipping'.format(model_file))
        return
    model_galaxy_mcmc(model_file, output_name=output_name, **_mcparams)


group_per_node = len(indices_to_fit) // size
remainder = len(indices_to_fit) % size

# Distribute remainder across the first `remainder` ranks
if rank < remainder:
    min_in_rank = rank * (group_per_node + 1)
    max_in_rank = min_in_rank + group_per_node + 1
else:
    min_in_rank = rank * group_per_node + remainder
    max_in_rank = min_in_rank + group_per_node

galaxies_cycle = indices_to_fit[min_in_rank:max_in_rank]

print(f"Rank {rank} processing indices {min_in_rank} to {max_in_rank}")

for count, ind_current in enumerate(galaxies_cycle):
    print(f"Processing BH #...{count} out of {len(galaxies_cycle)} on rank # {rank}")
    if config.noiseless_psf:
        model_file = f'/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_model_files_noiseless_{filter_first_str}_wo_quasar_fits_files_6p5_exp_{exp_time}_FOV_{FOVs}_samp_{ss}/without_quasar_model_{ind_current}.py'
    else:
        model_file = f'/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/PSFMC_MODEL_FILES/all_model_files_{filter_first_str}_wo_quasar_fits_files_6p5_exp_{exp_time}_FOV_{FOVs}_samp_{ss}/without_quasar_model_{ind_current}.py'

    print(f"Running model file {model_file}")
    try:
        run_mcmc(model_file)
    except Exception as E:
        print(E)
        continue




# for r in range(comm.rank, size):
#     for i in range(num_galaxies_per_task):
#         t_start = MPI.Wtime()
#
#         ind_current = (num_galaxies_per_task * comm.rank) + i  # don't want to multiply by 0
#
#         print("Processing BH #...")
#         print(ind_current)
#         print("on rank # ", comm.rank)
#         model_file = f'/fred/oz183/sberger/paper_2_obs_bias/all_model_files/without_quasar_model_{ind_current}.py'  # 'mcmc_model_mock_JWST_SDSS_{}.py'.format(ext))
#         run_mcmc(model_file)
#         t_diff = MPI.Wtime() - t_start  ### Stop stopwatch ###
#         print("This iteration took, ", t_diff, " seconds.")
