import numpy as np
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import scipy.stats as stats

# Step 1: Generate a bimodal dataset, which is a common example of a complex distribution
# a single Gaussian cannot model well.
np.random.seed(0)
# Data from the first Gaussian component
data1 = np.random.normal(-4, 1.0, 300)
# Data from the second Gaussian component
data2 = np.random.normal(4, 1.5, 700)
# Combine them into one dataset
data = np.concatenate((data1, data2)).reshape(-1, 1)

# Step 2: Fit a Gaussian Mixture Model with K=2 components
gmm = GaussianMixture(n_components=2, random_state=0)
gmm.fit(data)

# Step 3: Print the parameters of the fitted GMM.
# A GMM's probability density function is a weighted sum of Gaussian densities.
# The "equation" is p(x) = sum_{k=1 to K} [weight_k * N(x | mean_k, covariance_k)]
print("Fitted Gaussian Mixture Model Parameters:")
print("="*40)
for i in range(gmm.n_components):
    weight = gmm.weights_[i]
    mean = gmm.means_[i][0]
    cov = gmm.covariances_[i][0][0]
    print(f"Component {i+1}:")
    print(f"  Weight (π_{i+1}): {weight:.4f}")
    print(f"  Mean   (μ_{i+1}): {mean:.4f}")
    print(f"  Variance (σ²_{i+1}): {cov:.4f}")
print("="*40)
print("\nThe final model is a weighted sum of these two Gaussian distributions.")
print("This demonstrates that a mixture model can identify the underlying bimodal structure")
print("of the data, which a single Gaussian would fail to capture.")
print("\nThis supports option A: a mixture of Gaussians can model more complex, multi-modal distributions.")

# Optional: Visualization to show the superior fit
plt.figure(figsize=(10, 6))
plt.hist(data, bins=30, density=True, alpha=0.6, label='Data Histogram')

# Plot the single Gaussian fit (for comparison)
mu, std = stats.norm.fit(data)
x = np.linspace(data.min(), data.max(), 100)
p_single_gaussian = stats.norm.pdf(x, mu, std)
plt.plot(x, p_single_gaussian, 'r--', lw=2, label='Single Gaussian Fit (Poor)')

# Plot the GMM fit
p_gmm = np.exp(gmm.score_samples(x.reshape(-1, 1)))
plt.plot(x, p_gmm, 'g-', lw=2, label='Gaussian Mixture Model Fit (Good)')

plt.title('Fit Comparison: Single Gaussian vs. GMM')
plt.xlabel('Value')
plt.ylabel('Density')
plt.legend()
# To prevent the plot from showing in a headless environment and stopping the script
# we save it instead of showing it.
plt.savefig("gmm_fit_comparison.png")
print("\nA plot 'gmm_fit_comparison.png' has been saved to the current directory to visualize the fits.")
