import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# --- 1. Generate Synthetic Multi-Modal Data ---
# This data represents a situation where a single Gaussian is inaccurate.
np.random.seed(0)
# Data from the first mode (e.g., a sub-population)
data1 = np.random.normal(loc=-4, scale=1.5, size=300)
# Data from the second mode (e.g., another sub-population)
data2 = np.random.normal(loc=5, scale=2, size=700)
# Combine into one dataset
data = np.concatenate((data1, data2)).reshape(-1, 1)

# --- 2. Fit a Single Gaussian Distribution ---
# This is the inaccurate model the user wants to replace.
mu_single, std_single = norm.fit(data)

# --- 3. Fit a Gaussian Mixture Model (GMM) ---
# This is the proposed solution. We choose K=2 components because we
# can see (or suspect) two underlying groups.
gmm = GaussianMixture(n_components=2, random_state=0)
gmm.fit(data)

# --- 4. Prepare for Plotting ---
# Create a range of x values to plot the PDFs
x_plot = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

# Get the PDF for the single Gaussian
pdf_single = norm.pdf(x_plot, mu_single, std_single)

# Get the PDF for the GMM. The GMM's PDF is a weighted sum of the PDFs of its components.
log_pdf_gmm = gmm.score_samples(x_plot)
pdf_gmm = np.exp(log_pdf_gmm)

# --- 5. Generate and Display the Plot ---
print("Generating plot to demonstrate the difference between a single Gaussian and a GMM fit...")
print("Please close the plot window to continue.")

plt.figure(figsize=(12, 7))
# Plot the data histogram
plt.hist(data, bins=50, density=True, alpha=0.6, color='skyblue', label='Real-world Data Histogram')
# Plot the single Gaussian fit
plt.plot(x_plot, pdf_single, 'r--', linewidth=2, label='Single Gaussian Fit (Inaccurate)')
# Plot the GMM fit
plt.plot(x_plot, pdf_gmm, 'g-', linewidth=2, label='Gaussian Mixture Fit (K=2, Accurate)')

plt.title('GMM vs. Single Gaussian for Multi-Modal Data')
plt.xlabel('Observation Value')
plt.ylabel('Probability Density')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()

# --- 6. Print the Conclusion ---
# Explaining why the GMM is a good idea, based on the visualization
print("\n--- Analysis ---")
print("The plot demonstrates that the data has two distinct peaks (it is bimodal).")
print(f"The single Gaussian model (mean={mu_single:.2f}, std={std_single:.2f}) fails to capture this structure and provides a poor fit.")
print("The Gaussian Mixture Model, however, uses a combination of two Gaussians to accurately model both peaks.")
print("\nThis confirms that using a mixture of Gaussians is a good idea because it can model more complex, multi-modal distributions that are common in real-world data.")
print("This aligns with the reasoning in option A.")
