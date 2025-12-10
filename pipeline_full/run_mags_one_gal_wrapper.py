import subprocess
import numpy as np
import sys

filter_type = sys.argv[1]
samps = sys.argv[2]

# magnitudes_range = np.arange(22.1, 25, 0.1)
magnitudes_range = [22.1]
np.save("noiseless_test.npy", np.asarray([0, 1, 2, 3, 4, 5, 6]))
for magnitude in magnitudes_range:
    process = subprocess.run(["python", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/wrapper.py",
                    "--run_through_psfmc", "True",
                    "--include_quasar", "True",
                    "--generate_images", "False",
                    "--filters", filter_type,
                    "--mags_AB", str(np.round(magnitude, 1)),
                    "--samps", str(samps),
                    "--exp_time", str(10000),
                    "--noiseless_psf", "True",
                    "--indices_to_fit", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/noiseless_test.npy",
                    "--MPI_tasks", "3"])
