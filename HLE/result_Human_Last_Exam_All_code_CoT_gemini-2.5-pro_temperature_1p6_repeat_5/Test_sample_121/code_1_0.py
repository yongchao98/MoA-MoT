import numpy as np
from sklearn.mixture import GaussianMixture

# --- Step 1: Create synthetic bimodal data ---
# This data represents observations from two distinct sub-populations or modes.
np.random.seed(42) # for reproducibility
data_cluster1 = np.random.normal(loc=-5, scale=1.5, size=300)
data_cluster2 = np.random.normal(loc=5, scale=2.0, size=700)
data = np.concatenate([data_cluster1, data_cluster2]).reshape(-1, 1)

# --- Step 2: Fit a Gaussian Mixture Model ---
# We hypothesize there are two underlying modes, so we set n_components=2.
gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(data)

# --- Step 3: Extract and print the parameters ---
# A single Gaussian would fail to capture the two peaks. The GMM, however,
# creates an emission density that is a sum of two Gaussians.
# The equation for the density is:
# P(x) = weight_1 * Gaussian(x | mean_1, std_1) + weight_2 * Gaussian(x | mean_2, std_2)

print("Gaussian Mixture Model has been successfully fitted.")
print("The model found the following components to describe the data:\n")

# Get the parameters
weights = gmm.weights_
means = gmm.means_
# Covariances are variance for 1D data. Get std dev by taking the square root.
stds = np.sqrt(gmm.covariances_)

# --- Step 4: Output each number in the final equation ---
for i in range(gmm.n_components):
    print(f"--- Component {i+1} ---")
    print(f"This component represents a fraction of the data defined by the weight.")
    print(f"Weight (w_{i+1}): {weights[i]:.4f}")
    print("\nThis component is a Gaussian distribution with the following parameters:")
    print(f"Mean (μ_{i+1}): {means[i][0]:.4f}")
    print(f"Standard Deviation (σ_{i+1}): {stds[i][0]:.4f}\n")

print("The final emission density is a weighted sum of these two Gaussian distributions.")
print("This allows the model to accurately represent the complex, two-peaked nature of the data.")