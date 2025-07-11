import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Demonstrates that a Gaussian Mixture Model (GMM) provides a better fit
    for multi-modal data compared to a single Gaussian distribution.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Generate multi-modal data by combining two different Gaussian distributions
    # This simulates real-world data with underlying sub-populations.
    data1 = np.random.normal(loc=-4, scale=1.5, size=300)
    data2 = np.random.normal(loc=5, scale=2.0, size=700)
    data = np.concatenate([data1, data2]).reshape(-1, 1)

    print("--- Data Generated ---")
    print(f"Total data points: {len(data)}")
    print("Data is composed of two clusters, one around -4 and another around 5.\n")

    # --- Case 1: Fit a Single Gaussian Distribution ---
    print("--- Model 1: Single Gaussian ---")
    mean_single = np.mean(data)
    std_single = np.std(data)

    # The "equation" for a single Gaussian N(x | mu, sigma^2)
    print("Fitted Parameters (The numbers in the equation):")
    print(f"  - Mean (mu): {mean_single:.4f}")
    print(f"  - Std. Dev. (sigma): {std_single:.4f}\n")

    # Calculate log-likelihood for the single Gaussian model
    # Log-likelihood is a measure of how well the model explains the data. Higher is better.
    log_likelihood_single = np.sum(norm.logpdf(data, loc=mean_single, scale=std_single))
    print(f"Total Log-Likelihood for Single Gaussian Model: {log_likelihood_single:.4f}\n\n")


    # --- Case 2: Fit a Gaussian Mixture Model (GMM) ---
    print("--- Model 2: Gaussian Mixture Model (2 Components) ---")
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    # The "equation" for a GMM is a weighted sum of Gaussians:
    # p(x) = w1*N(x|mu1, sig1) + w2*N(x|mu2, sig2)
    print("Fitted Parameters (The numbers in the equations):")
    for i in range(gmm.n_components):
        print(f"  Component {i+1}:")
        print(f"    - Weight (w{i+1}): {gmm.weights_[i]:.4f}")
        print(f"    - Mean (mu{i+1}): {gmm.means_[i, 0]:.4f}")
        # sklearn's covariance is variance for 1D data
        print(f"    - Std. Dev. (sigma{i+1}): {np.sqrt(gmm.covariances_[i, 0]):.4f}")
    print()

    # The score_samples method gives log P(x_i | model), we sum them for total log-likelihood
    log_likelihood_gmm = np.sum(gmm.score_samples(data))
    print(f"Total Log-Likelihood for GMM: {log_likelihood_gmm:.4f}\n\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    if log_likelihood_gmm > log_likelihood_single:
        print("The GMM has a significantly higher log-likelihood, indicating a much better fit.")
        print("This is because the GMM can model the two separate clusters (multi-modality) in the data,")
        print("whereas the single Gaussian model fails to capture this complex structure.")
    else:
        print("The single Gaussian model provided a better or equal fit (this is highly unlikely for this data).")


if __name__ == '__main__':
    demonstrate_gmm_fit()