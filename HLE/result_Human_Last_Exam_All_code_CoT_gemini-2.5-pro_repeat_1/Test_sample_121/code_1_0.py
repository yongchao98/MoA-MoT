import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# --- Step 1: Create synthetic real-world data ---
# This data will have two distinct clusters (bimodal), which is common in real datasets.
# A single Gaussian distribution will be a poor fit for this.
np.random.seed(42)
# Create the first cluster of data points
cluster1 = np.random.normal(loc=-4, scale=1.0, size=300)
# Create the second cluster of data points
cluster2 = np.random.normal(loc=5, scale=1.5, size=700)
# Combine the clusters into a single dataset
data = np.concatenate((cluster1, cluster2)).reshape(-1, 1)

# --- Step 2: Fit a single Gaussian distribution ---
# This is the baseline approach with an inaccurate emission density.
mean_single = np.mean(data)
std_single = np.std(data)

# --- Step 3: Fit a Gaussian Mixture Model (GMM) ---
# This is the proposed, more accurate approach. We choose K=2 components
# because we know our synthetic data has two clusters.
gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(data)
weights = gmm.weights_
means = gmm.means_
# In 1D, covariance is just the variance (std_dev^2)
stds = np.sqrt(gmm.covariances_)

# --- Step 4: Print the results and model parameters ---
# We will output the "equation" for both models to show how they differ.

print("="*60)
print("Demonstrating the benefit of a Gaussian Mixture emission density")
print("="*60)

print("\n--- Model 1: Single Gaussian Fit ---")
print("This model attempts to fit all data with one bell curve.")
print("The probability density function (PDF) is: N(x | μ, σ²)")
print("\nFinal Equation (Parameters):")
print(f"N(x | μ={mean_single:.2f}, σ²={std_single**2:.2f})")
print("\nObservation: This single model poorly represents the two underlying")
print("data clusters by averaging them into one wide distribution.")


print("\n\n--- Model 2: Gaussian Mixture Model (GMM) Fit with K=2 ---")
print("This model fits the data with a weighted sum of 2 bell curves.")
print("The PDF is: w₁*N(x | μ₁, σ₁²) + w₂*N(x | μ₂, σ₂²)")
print("\nFinal Equation (Parameters):")
# We sort the components by their mean for consistent output
sorted_indices = np.argsort(means.flatten())
w1, w2 = weights[sorted_indices]
m1, m2 = means[sorted_indices].flatten()
s1_sq, s2_sq = gmm.covariances_[sorted_indices].flatten()

print(f"Component 1: weight={w1:.2f}, N(x | μ₁={m1:.2f}, σ₁²={s1_sq:.2f})")
print(f"Component 2: weight={w2:.2f}, N(x | μ₂={m2:.2f}, σ₂²={s2_sq:.2f})")
print(f"\nResulting PDF: {w1:.2f}*N(x|{m1:.2f},{s1_sq:.2f}) + {w2:.2f}*N(x|{m2:.2f},{s2_sq:.2f})")
print("\nObservation: The GMM correctly identifies the two modes (-4 and 5)")
print("and models them as separate components, providing a much more")
print("accurate and descriptive model of the complex data.")
