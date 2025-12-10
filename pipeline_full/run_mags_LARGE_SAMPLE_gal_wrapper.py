import subprocess
import numpy as np
import sys

filter_type = sys.argv[1]
samps = sys.argv[2]
quasar_mag = sys.argv[3]
indices_file = sys.argv[4]
exp_time = sys.argv[5]
psf_mag = sys.argv[6]

print(f"filter_type: {filter_type}")
print(f"samps: {samps}")
print(f"quasar_mag: {quasar_mag}")
print(f"indices_file: {indices_file}")
print(f"psf_mag: {psf_mag}")

# magnitudes_range = np.arange(22.1, 25, 0.1)
magnitudes_range = [quasar_mag]

for magnitude in magnitudes_range:
    process = subprocess.run(["python", "/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/wrapper.py",
                    "--run_through_psfmc", "True",
                    "--include_quasar", "True",
                    "--generate_images", "True",
                    "--filters", filter_type,
                    "--mags_AB", str(magnitude),
                    "--samps", str(samps),
                    "--exp_time", str(exp_time),
                    "--noiseless_psf", "False",
                    "--MPI_tasks", "10",
                    "--indices_to_fit", f"/fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/indices_to_fit_direc/{indices_file}",
                    "--full", "False",
                    "--psf_mag", psf_mag])


