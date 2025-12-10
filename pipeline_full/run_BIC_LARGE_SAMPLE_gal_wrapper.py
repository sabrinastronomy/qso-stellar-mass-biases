import subprocess
import numpy as np
import sys

filter_type = sys.argv[1]
samps = sys.argv[2]

# magnitudes_range = np.arange(22.1, 25, 0.1)
magnitudes_range = [19.543]

for magnitude in magnitudes_range:
    process = subprocess.run(["python", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/wrapper.py",
                    "--run_through_psfmc", "True",
                    "--include_quasar", "True",
                    "--generate_images", "False",
                    "--BIC", "True",
                    "--filters", filter_type,
                    "--mags_AB", str(np.round(magnitude, 1)),
                    "--samps", str(samps),
                    "--exp_time", str(10000),
                    "--noiseless_psf", "False",
                    "--MPI_tasks", "10",
                    "--indices_to_fit", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/selected_indices/FULL_selected_indices_sample_23.5_24.5.npy",
                    "--full", "True"])
