#!/bin/bash
#SBATCH --job-name=run_all_mags_J0148+0600_F115W
#SBATCH --time=10:00:00
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -D /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/
#SBATCH -o /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/out_err_files/run_all_mags_J0148+0600_F115W.out
#SBATCH -e /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/out_err_files/run_all_mags_J0148+0600_F115W.err
#SBATCH --mail-user=sabrinaberger55@gmail.com
#SBATCH --mail-type=ALL

source /home/sberger/.bash_profile
conda activate /home/sberger/.conda/envs/z-quasar
python /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/run_mags_LARGE_SAMPLE_gal_wrapper.py JWST.NIRCAM.F115W 1 19.522 BY_WIDTH_J0148+0600_F115W_photo_cropped_selected_indices_sample_21.48_25.48.npy 3100 17
