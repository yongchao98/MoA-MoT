import numpy as np
from sklearn.mixture import GaussianMixture

# Step 1: Create a synthetic multi-modal dataset.
# This data simulates a real-world scenario where observations come from two distinct
# sub-populations, making the overall distribution non-Gaussian.
np.random.seed(0)
# Create the first sub-population (cluster)
data_cluster1 = np.random.normal(loc=-4.0, scale=1.0, size=300)
# Create the second sub-population (cluster)
data_cluster2 = np.random.normal(loc=5.0, scale=1.5, size=700)
# Combine them into a single dataset
data = np.concatenate((data_cluster1, data_cluster2)).reshape(-1, 1)

# Step 2: Construct the emission density from a mixture of 2 Gaussians.
# We fit a Gaussian Mixture Model (GMM) with K=2 components to the data.
gmm = GaussianMixture(n_components=2, random_state=0)
gmm.fit(data)

# Step 3: Extract and print the learned parameters of the mixture model.
# A GMM models the data's probability density as a weighted sum of Gaussian distributions.
# The general equation is: P(x) = weight_1 * Gaussian_1(x) + weight_2 * Gaussian_2(x)
print("Demonstration: Fitting a GMM to multi-modal data.")
print("The generated data has two peaks, which a single Gaussian cannot model effectively.")
print("The GMM, however, can learn the parameters of the two underlying distributions.\n")

# Get the learned parameters from the fitted model
means = gmm.means_.flatten()
variances = gmm.covariances_.flatten()
weights = gmm.weights_.flatten()

# Ensure components are printed in a consistent order for clarity
# Sort by mean
sort_indices = np.argsort(means)
means = means[sort_indices]
variances = variances[sort_indices]
weights = weights[sort_indices]

# Print the parameters for each component that form the final "equation"
print("The model discovered the following two Gaussian components:")
print("-" * 55)
# Output numbers for the first Gaussian in the equation
print(f"Component 1 (approximating the first data cluster):")
print(f"  - Weight (w1): {weights[0]:.4f}")
print(f"  - Mean (mu1):   {means[0]:.4f}")
print(f"  - Variance (sigma^2_1): {variances[0]:.4f}")
print("-" * 55)
# Output numbers for the second Gaussian in the equation
print(f"Component 2 (approximating the second data cluster):")
print(f"  - Weight (w2): {weights[1]:.4f}")
print(f"  - Mean (mu2):   {means[1]:.4f}")
print(f"  - Variance (sigma^2_2): {variances[1]:.4f}")
print("-" * 55)

# Step 4: Display the final equation using the learned parameters.
# N(x | mu, sigma^2) represents a Gaussian probability density function.
print("\nThe final emission density can be represented by the equation:")
print(f"P(x) = {weights[0]:.4f} * N(x | {means[0]:.4f}, {variances[0]:.4f}) + {weights[1]:.4f} * N(x | {means[1]:.4f}, {variances[1]:.4f})")
print("\nThis demonstrates the GMM's ability to model a complex, multi-modal distribution.")