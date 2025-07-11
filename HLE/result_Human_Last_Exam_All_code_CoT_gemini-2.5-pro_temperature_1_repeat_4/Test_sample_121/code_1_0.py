import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# 1. Generate bimodal data
np.random.seed(42)
# Samples from the first Gaussian component
data1 = np.random.normal(loc=-4, scale=1, size=300)
# Samples from the second Gaussian component
data2 = np.random.normal(loc=4, scale=1.5, size=700)
# Combine into one dataset
data = np.concatenate([data1, data2]).reshape(-1, 1)

# 2. Fit a single Gaussian
mu_single = np.mean(data)
std_single = np.std(data)

# 3. Fit a Gaussian Mixture Model with 2 components
gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(data)

# 4. Print the parameters of the fitted GMM ("the final equation")
print("Fitted Gaussian Mixture Model Parameters:")
for i in range(gmm.n_components):
    weight = gmm.weights_[i]
    mean = gmm.means_[i][0]
    variance = gmm.covariances_[i][0][0]
    print(f"\nComponent {i+1}:")
    print(f"  Weight (mixing coefficient): {weight:.4f}")
    print(f"  Mean: {mean:.4f}")
    print(f"  Variance: {variance:.4f}")

print("\n---")
print("For comparison, the single Gaussian fit has:")
print(f"  Mean: {mu_single:.4f}")
print(f"  Variance: {std_single**2:.4f}")
print("\nPlotting the results...")


# 5. Visualize the results
plt.figure(figsize=(12, 7))
# Plot the histogram of the data
plt.hist(data, bins=30, density=True, alpha=0.6, color='gray', label='Real Data Histogram')

# Create a range of x values for plotting the PDFs
x_plot = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

# Plot the single Gaussian PDF
pdf_single = norm.pdf(x_plot, mu_single, std_single)
plt.plot(x_plot, pdf_single, 'r--', linewidth=2, label='Single Gaussian Fit')

# Plot the GMM PDF
# The score_samples method returns the log-likelihood, so we exponentiate it
log_pdf_gmm = gmm.score_samples(x_plot)
pdf_gmm = np.exp(log_pdf_gmm)
plt.plot(x_plot, pdf_gmm, 'b-', linewidth=2, label='Gaussian Mixture Model Fit (K=2)')

plt.title('GMM vs. Single Gaussian for Bimodal Data')
plt.xlabel('Value')
plt.ylabel('Density')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()
