import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def generate_chemistry(n_specimens, n_peaks):
    """
    A Python implementation of the core logic from the R script.
    It generates one group of specimens.
    """
    # Generate peaks from normal distributions
    # Each peak 'm' has a mean 'm'
    peaks = np.zeros((n_specimens, n_peaks))
    for m in range(n_peaks):
        peaks[:, m] = np.random.normal(loc=m + 1, scale=1, size=n_specimens)

    # Generate a single baseline for the whole group
    baseline = np.random.uniform(0, 1, n_peaks)

    # Subtract the same baseline from all specimens
    # This is equivalent to the R code's apply loop for a single treatment group
    peaks_baselined = peaks - baseline

    # Normalize each specimen (row)
    # The [:, np.newaxis] is for broadcasting
    row_sums = peaks_baselined.sum(axis=1)
    # Avoid division by zero, although unlikely here
    row_sums[row_sums == 0] = 1 
    normalized_peaks = peaks_baselined / row_sums[:, np.newaxis]
    
    return normalized_peaks

# --- Simulate Control Group ---
n_specimens_control = 360
n_peaks_shared = 340
controls = generate_chemistry(n_specimens_control, n_peaks_shared)

# --- Simulate Water Stress Group ---
n_specimens_stress = 510
water_stress = generate_chemistry(n_specimens_stress, n_peaks_shared)

# --- Perform PCA and Visualize ---
pca_control = PCA(n_components=2)
controls_pca = pca_control.fit_transform(controls)

pca_stress = PCA(n_components=2)
stress_pca = pca_stress.fit_transform(water_stress)

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.scatter(controls_pca[:, 0], controls_pca[:, 1], alpha=0.6)
ax1.set_title("Control Group PCA")
ax1.set_xlabel("PC1")
ax1.set_ylabel("PC2")
ax1.grid(True)

ax2.scatter(stress_pca[:, 0], stress_pca[:, 1], alpha=0.6, color='orange')
ax2.set_title("Water Stress Group PCA")
ax2.set_xlabel("PC1")
ax2.set_ylabel("PC2")
ax2.grid(True)

plt.suptitle("Simulated PCA Plots")
plt.tight_layout()
# The plot will show as output when the code is executed.
# This visualization confirms that each group forms a single cluster.
plt.show()

# The final answer is based on the logical analysis of the R script.
# The control group has 1 chemotype.
# The water stressed group has 1 chemotype.
print("Based on the analysis of the R script's generative process:")
print("Expected number of chemotypes for the control group: 1")
print("Expected number of chemotypes for the water stressed group: 1")
