#!/bin/bash

# Load the Python dictionary from the file
python_dict_file="quasar_dict.py"

# Extract the dictionary contents into variables
readarray -t lines < <(python -c "
import sys
sys.path.insert(0, '.')
namespace = {}
with open('quasar_dict.py') as f:
    exec(f.read(), namespace)

def get_wavelength(filter_name):
    wavelengths = {
        'F115W': 1.15,
        'F150W': 1.50,
        'F200W': 2.00,
        'F356W': 3.56
    }
    return wavelengths.get(filter_name, 0)

quasar_data = namespace['quasar_data']
exposure_times_dict = namespace['exposure_times_dict']

for quasar, types in quasar_data.items():
    for obj_type, filters in types.items():
        is_galaxy = obj_type == 'galaxy'
        for filter, data in (filters or {}).items():
            if data:
                mag_value = float(data[0]) if is_galaxy else data[0]
                exp_time = exposure_times_dict.get('EIGER', {}).get(filter, exposure_times_dict.get('Ding', {}).get(filter, 0))
                wavelength = get_wavelength(filter)
                ss_value = 2 if wavelength > 3 else 1
                galaxy_mag_min = round(mag_value - 2, 2) if is_galaxy else None
                galaxy_mag_max = round(mag_value + 2, 2) if is_galaxy else None
                print(f'{quasar} {filter} {mag_value} {exp_time} {ss_value} {galaxy_mag_min} {galaxy_mag_max}')
")

# Create SLURM scripts for each combination
output_dir="slurm_scripts"
mkdir -p "$output_dir"

for line in "${lines[@]}"; do
    # Parse the data
    read -r quasar filter mag_value exp_time ss galaxy_mag_min galaxy_mag_max <<< "$line"

    # Format the quasar name and filter for file naming
    quasar_name="${quasar%%[+-]*}"
    quasar_filter="${filter}"

    # Create the SLURM script
    script_path="$output_dir/run_all_just_gal_${quasar_name}_${quasar_filter}.sh"

    cat <<EOL > "$script_path"
#!/bin/bash
#SBATCH --job-name=run_all_just_gal_${quasar_name}_${quasar_filter}
#SBATCH --time=10:00:00
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=70G
#SBATCH -D /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/
#SBATCH -o /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/slurm_scripts/out_err_files/run_all_just_gal_${quasar_name}_${quasar_filter}.out
#SBATCH -e /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/slurm_scripts/out_err_files/run_all_just_gal_${quasar_name}_${quasar_filter}.err
#SBATCH --mail-user=sabrinaberger55@gmail.com
#SBATCH --mail-type=ALL

source /home/sberger/.bash_profile
conda activate /home/sberger/.conda/envs/z-quasar
python /fred/oz183/sberger/paper_2_obs_bias/src/pipeline_full/run_just_gal_LARGE_SAMPLE_gal_wrapper.py JWST.NIRCAM.${filter} ${ss} BY_WIDTH_${quasar_name}_${filter}_all_indices_sample_${galaxy_mag_min}_${galaxy_mag_max}.npy ${exp_time} 17
EOL

    echo "Generated: $script_path"
done
