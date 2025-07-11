import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# 1. Generate a synthetic multi-modal dataset
# This represents a real-world scenario where data comes from more than one source/mode.
np.random.seed(0)
# Data from the first mode
data1 = np.random.normal(-3, 1, 300)
# Data from the second mode
data2 = np.random.normal(4, 1.5, 700)
# Combine them into one dataset
data = np.concatenate((data1, data2)).reshape(-1, 1)

# 2. Fit a single Gaussian to the data
# This is the less accurate model that we want to improve upon.
mean_single = data.mean()
std_single = data.std()
# The "equation" for this single Gaussian is N(mean, std^2)
print("--- Single Gaussian Fit ---")
print(f"This model assumes all data comes from a single distribution:")
print(f"Mean = {mean_single:.4f}")
print(f"Standard Deviation = {std_single:.4f}")
print("-" * 29 + "\n")


# 3. Fit a Gaussian Mixture Model (GMM) with K=2 components
# This is the proposed, more accurate model.
gmm = GaussianMixture(n_components=2, random_state=0)
gmm.fit(data)

# Extract GMM parameters to show the "equation" for the mixture
weights = gmm.weights_
means = gmm.means_
covariances = gmm.covariances_

print("--- Gaussian Mixture Model Fit (K=2) ---")
print("This model finds a mixture of two distinct distributions:")
for i in range(len(weights)):
    print(f"Component {i+1}:")
    # Each number in the final equation is printed below
    print(f"  Weight (proportion of data): {weights[i]:.4f}")
    print(f"  Mean: {means[i][0]:.4f}")
    print(f"  Standard Deviation: {np.sqrt(covariances[i][0][0]):.4f}")
print("-" * 42 + "\n")


# 4. Visualize the results to show GMM is a better fit
plt.figure(figsize=(12, 6))
# Plot the histogram of the data
plt.hist(data, bins=50, density=True, alpha=0.6, color='gray', label='Observed Data')

# Create a range of x-values for plotting the PDFs
x_axis = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

# Plot the PDF of the single Gaussian fit
pdf_single = norm.pdf(x_axis, mean_single, std_single)
plt.plot(x_axis, pdf_single, 'r--', linewidth=2, label='Single Gaussian Fit')

# Plot the PDF of the GMM fit
# The total GMM PDF is the weighted sum of individual component PDFs
pdf_gmm = np.exp(gmm.score_samples(x_axis))
plt.plot(x_axis, pdf_gmm, 'b-', linewidth=2, label='Gaussian Mixture Model Fit')

plt.title('GMM vs. Single Gaussian for Modeling Multi-Modal Data')
plt.xlabel('Value')
plt.ylabel('Density')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

print("The plot demonstrates that the GMM (blue line) accurately captures the two modes of the data,")
print("while the single Gaussian (red dashed line) fails to do so, providing a poor overall fit.")
print("This visually confirms the principle described in option A.")

plt.show()
