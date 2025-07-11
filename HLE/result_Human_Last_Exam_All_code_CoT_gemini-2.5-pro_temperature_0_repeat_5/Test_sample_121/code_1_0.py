import numpy as np
from sklearn.mixture import GaussianMixture

# 1. Generate a bimodal dataset
# Let's create data from two distinct Gaussian distributions to simulate
# a complex, multi-modal emission distribution.
np.random.seed(0)
# Component 1: 300 samples from a Gaussian with mean -4 and standard deviation 1
data1 = np.random.normal(loc=-4, scale=1, size=300)
# Component 2: 300 samples from a Gaussian with mean 4 and standard deviation 1.5
data2 = np.random.normal(loc=4, scale=1.5, size=300)

# Combine into a single dataset
data = np.concatenate([data1, data2]).reshape(-1, 1)

# 2a. Fit a single Gaussian model
gmm_single = GaussianMixture(n_components=1, random_state=0)
gmm_single.fit(data)

# 2b. Fit a Gaussian Mixture Model with 2 components
gmm_mixture = GaussianMixture(n_components=2, random_state=0)
gmm_mixture.fit(data)

# 3. Print and compare the results
print("--- Model 1: Single Gaussian ---")
print(f"This model assumes all data comes from one source.")
mean_single = gmm_single.means_[0][0]
cov_single = gmm_single.covariances_[0][0][0]
# The 'equation' for this density is N(mu, sigma^2)
print(f"Fitted Parameters (Equation): N(mu={mean_single:.4f}, sigma^2={cov_single:.4f})")
# The score() method gives average log-likelihood. Multiply by N for total.
log_likelihood_single = gmm_single.score(data) * data.shape[0]
print(f"Total Log-Likelihood: {log_likelihood_single:.4f}\n")


print("--- Model 2: Gaussian Mixture Model (2 Components) ---")
print("This model allows data to come from a mixture of sources.")
# The 'equation' for this density is w1*N(mu1, sigma1^2) + w2*N(mu2, sigma2^2)
weights = gmm_mixture.weights_
means = gmm_mixture.means_
covs = gmm_mixture.covariances_

print("Fitted Parameters (Equation):")
for i in range(len(weights)):
    print(f"  Component {i+1}: Weight(w{i+1})={weights[i]:.4f}, N(mu{i+1}={means[i][0]:.4f}, sigma{i+1}^2={covs[i][0][0]:.4f})")

log_likelihood_mixture = gmm_mixture.score(data) * data.shape[0]
print(f"Total Log-Likelihood: {log_likelihood_mixture:.4f}\n")

print("--- Conclusion ---")
print("The Gaussian Mixture Model has a significantly higher log-likelihood, indicating a much better fit.")
print("Its parameters also successfully recovered the two original distributions we used to generate the data.")
print("This demonstrates that using a mixture of Gaussians is a good idea to model complex, multi-modal data.")
