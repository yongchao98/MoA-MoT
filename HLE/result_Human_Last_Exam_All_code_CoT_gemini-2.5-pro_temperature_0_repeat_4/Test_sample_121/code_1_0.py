import numpy as np
from sklearn.mixture import GaussianMixture

def demonstrate_gmm_fit():
    """
    Demonstrates that a Gaussian Mixture Model provides a better fit
    to multi-modal data than a single Gaussian distribution.
    """
    # 1. Generate synthetic multi-modal data
    # Data from a mixture of two Gaussians
    np.random.seed(0)
    data1 = np.random.normal(loc=-4, scale=1, size=300)
    data2 = np.random.normal(loc=4, scale=1.5, size=700)
    data = np.concatenate([data1, data2]).reshape(-1, 1)

    # 2. Fit a single Gaussian model (GMM with 1 component)
    single_gaussian = GaussianMixture(n_components=1, random_state=0)
    single_gaussian.fit(data)
    # The score is the per-sample average log-likelihood.
    log_likelihood_single = single_gaussian.score(data)

    # 3. Fit a Gaussian Mixture Model (with 2 components)
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)
    log_likelihood_gmm = gmm.score(data)

    # 4. Print the results
    print("Demonstrating model fit on multi-modal data:")
    print(f"Log-Likelihood of Single Gaussian Model: {log_likelihood_single:.4f}")
    print(f"Log-Likelihood of 2-Component GMM: {log_likelihood_gmm:.4f}")
    print("\nConclusion: The GMM has a significantly higher log-likelihood, indicating a much better fit for the data.")
    print("This supports the idea that using a mixture of Gaussians is beneficial for modeling complex, non-Gaussian data.")

if __name__ == '__main__':
    demonstrate_gmm_fit()