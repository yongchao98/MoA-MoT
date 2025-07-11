import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# 1. Generate a synthetic bimodal dataset
# This represents data from a single state that doesn't follow a simple Gaussian distribution.
np.random.seed(42)
data_part1 = np.random.normal(loc=-2, scale=0.5, size=300)
data_part2 = np.random.normal(loc=3, scale=1.0, size=700)
data = np.concatenate([data_part1, data_part2]).reshape(-1, 1)

# 2. Fit a single Gaussian model
mean_sg = np.mean(data)
std_sg = np.std(data)
# Calculate log-likelihood for the single Gaussian model
log_likelihood_sg = np.sum(norm.logpdf(data, loc=mean_sg, scale=std_sg))

print("--- Single Gaussian Fit ---")
print("This model assumes all data comes from one simple distribution.")
# Equation for a single Gaussian: N(x | mu, sigma^2)
print(f"Fitted Mean (mu): {mean_sg:.4f}")
print(f"Fitted Std Dev (sigma): {std_sg:.4f}")
print(f"Total Log-Likelihood: {log_likelihood_sg:.4f}\n")


# 3. Fit a Gaussian Mixture Model (GMM) with 2 components
gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(data)
log_likelihood_gmm = gmm.score(data) * len(data) # .score() returns average log-likelihood

print("--- Gaussian Mixture Model (K=2) Fit ---")
print("This model allows for multiple underlying distributions (sub-populations).")
# Equation for a GMM: sum_{k=1 to K} [pi_k * N(x | mu_k, sigma_k^2)]

# Get and sort parameters for consistent output
weights = gmm.weights_
means = gmm.means_
covs = gmm.covariances_

# Sort by means to have a consistent order
sorted_indices = np.argsort(means.flatten())
weights = weights[sorted_indices]
means = means[sorted_indices]
covs = covs[sorted_indices]

# Output the parameters for each component in the final "equation"
for i in range(gmm.n_components):
    print(f"\nComponent {i+1}:")
    print(f"Weight (pi_{i+1}): {weights[i]:.4f}")
    print(f"Mean (mu_{i+1}): {means[i][0]:.4f}")
    print(f"Std Dev (sigma_{i+1}): {np.sqrt(covs[i][0][0]):.4f}")

print(f"\nTotal Log-Likelihood: {log_likelihood_gmm:.4f}")
print("\nConclusion: The GMM has a significantly higher log-likelihood, indicating a much better fit for this complex data.")
