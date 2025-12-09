#!/bin/bash
#SBATCH --job-name=run_all_mags_J0844-0132_F356W
#SBATCH --time=10:00:00
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -D /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/
#SBATCH -o /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/out_err_files/run_all_mags_J0844-0132_F356W.out
#SBATCH -e /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/out_err_files/run_all_mags_J0844-0132_F356W.err
#SBATCH --mail-user=sabrinaberger55@gmail.com
#SBATCH --mail-type=ALL

source /home/sberger/.bash_profile
conda activate /home/sberger/.conda/envs/z-quasar
python /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/run_mags_LARGE_SAMPLE_gal_wrapper.py JWST.NIRCAM.F356W 2 22.13 BY_WIDTH_J0844-0132_F356W_photo_cropped_selected_indices_sample_23.7_27.7.npy 3100 17
