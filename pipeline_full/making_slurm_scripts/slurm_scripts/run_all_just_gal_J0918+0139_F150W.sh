#!/bin/bash
#SBATCH --job-name=run_all_just_gal_J0918+0139_F150W
#SBATCH --time=10:00:00
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH -D /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/
#SBATCH -o /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/out_err_files/run_all_just_gal_J0918+0139_F150W.out
#SBATCH -e /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/out_err_files/run_all_just_gal_J0918+0139_F150W.err
#SBATCH --mail-user=sabrinaberger55@gmail.com
#SBATCH --mail-type=ALL

source /home/sberger/.bash_profile
conda activate /home/sberger/.conda/envs/z-quasar
python /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/run_just_gal_LARGE_SAMPLE_gal_wrapper.py JWST.NIRCAM.F150W 1 BY_WIDTH_J0918+0139_F150W_photo_cropped_selected_indices_sample_24.03_28.03.npy 3100 17
