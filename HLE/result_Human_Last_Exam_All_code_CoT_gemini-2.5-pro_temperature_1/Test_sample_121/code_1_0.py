import numpy as np
from sklearn.mixture import GaussianMixture

def demonstrate_gmm_fit():
    """
    This function demonstrates that a Gaussian Mixture Model (GMM) can fit complex,
    multi-modal data better than a single Gaussian model.
    """
    # 1. Generate synthetic multi-modal data
    # We create data from two different Gaussian distributions and combine them.
    # This simulates a real-world scenario where data has multiple sub-populations.
    np.random.seed(0)
    # Component 1: 300 samples, mean=0, std=1
    data1 = np.random.normal(0, 1, 300)
    # Component 2: 300 samples, mean=5, std=1.5
    data2 = np.random.normal(5, 1.5, 300)
    # Combine into one dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    print("Generated a multi-modal dataset with two underlying distributions.")
    print("-" * 50)

    # 2. Fit a single Gaussian model
    # This is equivalent to a GMM with n_components=1
    single_gaussian = GaussianMixture(n_components=1, random_state=0)
    single_gaussian.fit(data)
    # The score method returns the average log-likelihood of the data under the model.
    # A higher (less negative) score indicates a better fit.
    log_likelihood_single = single_gaussian.score(data)

    print("Fitting a single Gaussian model...")
    print(f"Mean of the single Gaussian: {single_gaussian.means_[0][0]:.4f}")
    print(f"Variance of the single Gaussian: {single_gaussian.covariances_[0][0][0]:.4f}")
    print(f"Average Log-Likelihood (Single Gaussian): {log_likelihood_single:.4f}")
    print("-" * 50)

    # 3. Fit a Gaussian Mixture Model with K=2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)
    log_likelihood_gmm = gmm.score(data)

    print("Fitting a Gaussian Mixture Model (GMM) with 2 components...")
    print(f"Means of the GMM components: {gmm.means_[0][0]:.4f}, {gmm.means_[1][0]:.4f}")
    print(f"Variances of the GMM components: {gmm.covariances_[0][0][0]:.4f}, {gmm.covariances_[1][0][0]:.4f}")
    print(f"Average Log-Likelihood (GMM): {log_likelihood_gmm:.4f}")
    print("-" * 50)

    # 4. Conclusion
    if log_likelihood_gmm > log_likelihood_single:
        print("Conclusion: The GMM provides a significantly better fit to the multi-modal data.")
        print("This supports the idea that using a mixture of Gaussians is beneficial for complex distributions.")
    else:
        print("Conclusion: The single Gaussian model performed better, which is unexpected.")

if __name__ == '__main__':
    demonstrate_gmm_fit()