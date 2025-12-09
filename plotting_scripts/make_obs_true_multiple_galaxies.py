"""
Make a truth vs observed galaxy magnitude plot for different galaxies and quasar magnitudes
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

truth_mag = np.load("truth_mag.npy")
print(truth_mag)
# for mag in ["22p5", "22p2", "22p7"]:
for mag in ["18p0", "21p0", "22p6", "22p7", "22p78", "22p8", "22p9", "24p0"]:
    observed_mag = np.load(f"observed_mag_{mag}.npy")
    ind = int(min(len(observed_mag), len(truth_mag)))
    print(ind)
    # ind = 2
    # res = linregress(observed_mag, truth_mag)
    truth_mag_sigs = np.load("truth_mag_sigs.npy")[:ind]
    observed_mag_sigs = np.load(f"observed_mag_sigs_{mag}.npy")[:ind]
    print(truth_mag[:ind])
    print(observed_mag[:ind])
    plt.scatter(truth_mag[:ind], truth_mag[:ind]-observed_mag[:ind], label=f"quasar mag = {mag}")
    # print(truth_mag[:ind])
    # plt.errorbar(observed_mag[:ind], truth_mag[:ind], xerr=observed_mag_sigs, alpha=0.8, label=f"quasar mag = {mag}", fmt="o")
    # plt.plot(observed_mag, res.intercept + res.slope*observed_mag, color="k", label="fit")

# plt.scatter(truth_mag, truth_mag, label="truth gal mag", c="black")
# for mag in ["21p0", "22p6", "22p7", "22p78", "22p8", "22p9", "24p0"]:

# quasar_mag = mag.replace("p", ".")

# ["21p0", "22p6", "22p7", "22p78", "22p8", "22p9", "24p0"]


plt.legend()
# plt.errorbar()
plt.xlabel(r"$\rm gal~m_{F150W}$ (TRUTH)")
plt.ylabel("residual (truth-observed)")
# plt.ylabel(r"$\rm gal~m_{F150W}$ (truth)")
plt.savefig("sample_comparison_mag.png", dpi=200)