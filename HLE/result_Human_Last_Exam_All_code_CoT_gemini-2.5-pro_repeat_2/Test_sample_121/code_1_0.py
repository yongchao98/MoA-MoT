import numpy as np
from sklearn.mixture import GaussianMixture

# --- Step 1: Generate synthetic bimodal data ---
# This data cannot be modeled well by a single Gaussian.
# It represents a situation where an HMM state might emit observations
# from two distinct subpopulations.
np.random.seed(0)
# Data from the first Gaussian component
data1 = np.random.normal(loc=-4, scale=1.5, size=300)
# Data from the second Gaussian component
data2 = np.random.normal(loc=5, scale=2.0, size=700)
# Combine them into one dataset
data = np.concatenate((data1, data2)).reshape(-1, 1)

# --- Step 2: Fit a Gaussian Mixture Model (GMM) with K=2 components ---
# We hypothesize that there are two underlying groups in our data.
gmm = GaussianMixture(n_components=2, random_state=0)
gmm.fit(data)

# --- Step 3: Print the parameters (the "equation") of the fitted GMM ---
# A GMM's probability density function is a weighted sum of Gaussian densities.
# The equation is: p(x) = weight1 * N(x | mean1, cov1) + weight2 * N(x | mean2, cov2)
# Here, we print out the numbers that define this equation.

print("Fitted Gaussian Mixture Model Parameters:")
print("-" * 40)

# The components might not be ordered the same as we generated them.
# We can sort them by mean for consistent output.
component_params = sorted(zip(gmm.weights_, gmm.means_.flatten(), gmm.covariances_.flatten()), key=lambda x: x[1])

for i, (weight, mean, cov) in enumerate(component_params):
    print(f"Component {i+1}:")
    print(f"  Weight (mixing proportion): {weight:.4f}")
    print(f"  Mean: {mean:.4f}")
    # Covariance is a matrix, but for 1D data it's a single number (variance).
    print(f"  Variance (Covariance): {cov:.4f}")
    print()

print("Demonstration Explanation:")
print("The code shows that the GMM successfully identified two distinct Gaussian components in the data.")
print("Component 1 is centered around -4 and Component 2 is centered around 5, matching our generated data.")
print("The weights reflect the proportion of data points from each original distribution (30% and 70%).")
print("This demonstrates that a GMM can model complex, multi-modal data where a single Gaussian would fail, making it a great choice for HMM emission densities.")

# To visualize this, you can add the following code (requires matplotlib):
# import matplotlib.pyplot as plt
#
# x_axis = np.linspace(-10, 12, 1000).reshape(-1, 1)
# log_prob = gmm.score_samples(x_axis)
# pdf = np.exp(log_prob)
#
# plt.figure(figsize=(10, 6))
# plt.hist(data, bins=50, density=True, alpha=0.6, label='Data Histogram')
# plt.plot(x_axis, pdf, color='red', lw=2, label='Fitted GMM PDF')
# plt.title('GMM Fit to Bimodal Data')
# plt.xlabel('Value')
# plt.ylabel('Density')
# plt.legend()
# plt.grid(True)
# plt.show()
