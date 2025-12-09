import numpy as np

# Load data from updated paths
# observed_mag = np.load("17psf_collated_data/observed_mag_test_356W.npy")
# observed_mag_likelihoods = np.load("17psf_collated_data/observed_mag_likelihoods_test_356W.npy")
# observed_mag_sigs = np.load("17psf_collated_data/observed_mag_sigs_test_356W.npy")
#
# truth_mag_wq = np.load("17psf_collated_data/truth_mag_wq_test_356W.npy")
# truth_mag_wq_likelihoods = np.load("17psf_collated_data/truth_mag_wq_likelihoods_test_356W.npy")
# truth_mag_sigs_wq = np.load("17psf_collated_data/truth_mag_sigs_wq_test_356W.npy")
#
# indices_used = np.load("17psf_collated_data/indices_used_test_356W.npy")

observed_mag = np.load("17psf_collated_data/observed_mag_J159-02_356W.npy")
observed_mag_likelihoods = np.load("17psf_collated_data/observed_mag_likelihoods_J159-02_356W.npy")
observed_mag_sigs = np.load("17psf_collated_data/observed_mag_sigs_J159-02_356W.npy")

truth_mag_wq = np.load("17psf_collated_data/truth_mag_wq_J159-02_356W.npy")
truth_mag_wq_likelihoods = np.load("17psf_collated_data/truth_mag_wq_likelihoods_J159-02_356W.npy")
truth_mag_sigs_wq = np.load("17psf_collated_data/truth_mag_sigs_wq_J159-02_356W.npy")

indices_used = np.load("17psf_collated_data/indices_used_J159-02_356W.npy")


#### getting smallest residual
print(f"Number of elements in truth_mag_wq: {len(truth_mag_wq)}")

# Compute absolute difference
diff = np.abs(truth_mag_wq - observed_mag)

# Find index of minimum difference
min_index = np.argmin(diff)

# Print results neatly
print(f"Minimum absolute difference: {diff[min_index]:.4f}")
print(f"Corresponding index used: {indices_used[min_index]}")
print(f"Truth magnitude: {truth_mag_wq[min_index]:.4f}")
print(f"Observed magnitude: {observed_mag[min_index]:.4f}")


print(f"Length of truth_mag_wq: {len(truth_mag_wq)}")
print(f"Length of observed_mag: {len(observed_mag)}")
print(f"Length of diff: {len(diff)}")
print(f"Length of indices_used: {len(indices_used)}")

print(f"Number of elements in truth_mag_wq: {len(truth_mag_wq)}")



# Get indices that would sort the array
sorted_indices = np.argsort(diff)

# Get the index of the second smallest difference
second_min_index = sorted_indices[1]

#### getting smallest residual
print(f"Number of elements in truth_mag_wq: {len(truth_mag_wq)}")

# Get indices that would sort the array
sorted_indices = np.argsort(diff)

# Get the index of the second smallest difference
second_min_index = sorted_indices[1]

# Print results neatly
print(f"Second absolute difference: {diff[second_min_index]:.4f}")
print(f"Corresponding index used: {indices_used[second_min_index]}")
print(f"Truth magnitude: {truth_mag_wq[second_min_index]:.4f}")
print(f"Observed magnitude: {observed_mag[second_min_index]:.4f}")

