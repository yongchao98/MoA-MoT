import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Demonstrates that a Gaussian Mixture Model (GMM) provides a better fit
    for complex, multi-modal data compared to a single Gaussian distribution.
    """
    # 1. Generate synthetic bimodal data
    np.random.seed(42) # for reproducibility
    # Data from first mode
    data1 = np.random.normal(loc=-4, scale=1, size=300)
    # Data from second mode
    data2 = np.random.normal(loc=4, scale=1.5, size=700)
    # Combine into a single dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    print("--- Data Generation ---")
    print(f"Generated a bimodal dataset with {len(data)} points.\n")

    # 2. Fit a single Gaussian model
    print("--- Model 1: Single Gaussian ---")
    # Fit the model (calculate mean and standard deviation)
    mu, std = norm.fit(data)
    # Calculate the total log-likelihood
    log_likelihood_single_gaussian = np.sum(norm.logpdf(data, loc=mu, scale=std))
    print(f"Fitted a single Gaussian with mean={mu:.2f} and std={std:.2f}")
    print(f"Total Log-Likelihood: {log_likelihood_single_gaussian:.2f}\n")

    # 3. Fit a Gaussian Mixture Model (GMM) with 2 components
    print("--- Model 2: Gaussian Mixture Model (K=2) ---")
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)
    # The 'score' method returns the average log-likelihood per sample
    # Multiply by the number of samples to get the total log-likelihood
    log_likelihood_gmm = gmm.score(data) * len(data)
    print("Fitted a GMM with 2 components.")
    for i in range(gmm.n_components):
        print(f"  Component {i+1}: Weight={gmm.weights_[i]:.2f}, Mean={gmm.means_[i,0]:.2f}, StdDev={np.sqrt(gmm.covariances_[i,0,0]):.2f}")
    print(f"Total Log-Likelihood: {log_likelihood_gmm:.2f}\n")

    # 4. Conclusion
    print("--- Conclusion ---")
    print("Comparing the total log-likelihoods:")
    print(f"Single Gaussian: {log_likelihood_single_gaussian:.2f}")
    print(f"GMM (K=2):      {log_likelihood_gmm:.2f}")

    if log_likelihood_gmm > log_likelihood_single_gaussian:
        print("\nThe GMM has a significantly higher log-likelihood, indicating a much better fit.")
        print("This is because the GMM can model the two underlying modes in the data,")
        print("whereas the single Gaussian cannot. This demonstrates the key advantage")
        print("described in option A: the ability to model complex, multi-modal distributions.")
    else:
        print("\nThe single Gaussian provided a better fit, which is unexpected for this data.")


if __name__ == '__main__':
    demonstrate_gmm_fit()