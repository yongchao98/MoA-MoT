import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# 1. Generate multi-modal data
# This data has two distinct peaks (bimodal), which is common in real-world scenarios.
np.random.seed(0)
# First component/cluster centered at -4
data1 = np.random.normal(-4, 0.8, 300)
# Second component/cluster centered at 3
data2 = np.random.normal(3, 1.2, 700)
data = np.concatenate((data1, data2)).reshape(-1, 1)

# 2. Fit a single Gaussian distribution
# This attempts to model the complex data with a single bell curve.
mu, std = norm.fit(data)
print("--- Single Gaussian Fit ---")
print(f"This model assumes all data comes from one source.")
print(f"Mean (μ): {mu:.4f}")
print(f"Standard Deviation (σ): {std:.4f}")
# The final equation for the single Gaussian probability density function (PDF) is:
# P(x) = (1 / (std * sqrt(2π))) * exp(-0.5 * ((x - mu) / std)²)
print(f"Final Equation Parameters: mu = {mu:.4f}, std = {std:.4f}\n")


# 3. Fit a Gaussian Mixture Model (GMM) with K=2 components
# This allows the model to find two underlying distributions.
gmm = GaussianMixture(n_components=2, random_state=0)
gmm.fit(data)

# Extract GMM parameters
weights = gmm.weights_
means = gmm.means_
covariances = gmm.covariances_

print("--- Gaussian Mixture Model Fit (K=2) ---")
print("This model assumes the data is a mix of two different sources.")
print("Component 1:")
print(f"  Weight: {weights[0]:.4f}")
print(f"  Mean (μ₁): {means[0][0]:.4f}")
print(f"  Std. Dev. (σ₁): {np.sqrt(covariances[0][0][0]):.4f}")
print("\nComponent 2:")
print(f"  Weight: {weights[1]:.4f}")
print(f"  Mean (μ₂): {means[1][0]:.4f}")
print(f"  Std. Dev. (σ₂): {np.sqrt(covariances[1][0][0]):.4f}")

# The final equation for the GMM PDF is a sum of the individual component PDFs:
# P(x) = weight₁ * PDF(x; μ₁, σ₁) + weight₂ * PDF(x; μ₂, σ₂)
print("\nFinal Equation Parameters:")
print(f"weight1 = {weights[0]:.4f}, mu1 = {means[0][0]:.4f}, std1 = {np.sqrt(covariances[0][0][0]):.4f}")
print(f"weight2 = {weights[1]:.4f}, mu2 = {means[1][0]:.4f}, std2 = {np.sqrt(covariances[1][0][0]):.4f}")

# 4. (Optional) Visualization code for the user to run
# This part of the code will generate a plot if run by the user,
# visually confirming that the GMM is a much better fit.
print("\n[INFO] Run this code in an environment with matplotlib to see the plot.")
try:
    x_axis = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

    # PDF for the single Gaussian
    pdf_single_gaussian = norm.pdf(x_axis, mu, std)

    # PDF for the GMM
    log_pdf_gmm = gmm.score_samples(x_axis)
    pdf_gmm = np.exp(log_pdf_gmm)

    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=30, density=True, alpha=0.6, color='gray', label='Data Histogram')
    plt.plot(x_axis, pdf_single_gaussian, 'r-', lw=2, label='Single Gaussian Fit')
    plt.plot(x_axis, pdf_gmm, 'g--', lw=3, label='Gaussian Mixture Model Fit (K=2)')
    plt.title('Single Gaussian vs. GMM Fit for Bimodal Data')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    # plt.show() # Uncomment to display the plot
except ImportError:
    print("[WARN] matplotlib is not installed. Skipping plot generation.")
