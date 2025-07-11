import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Demonstrates that a Gaussian Mixture Model (GMM) can model complex,
    multi-modal data far better than a single Gaussian distribution.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Create a bimodal dataset by combining two different Gaussian distributions
    # This simulates real-world data that doesn't fit a single bell curve.
    data1 = np.random.normal(loc=-4, scale=1, size=300)
    data2 = np.random.normal(loc=4, scale=1.5, size=700)
    data = np.concatenate([data1, data2]).reshape(-1, 1)
    
    print("--- Data Generation ---")
    print(f"Generated a dataset with {len(data)} points from two different distributions.")
    print("This creates a 'bimodal' or 'two-peaked' dataset.\n")

    # 2. Fit a single Gaussian model to the data
    mean = data.mean()
    std = data.std()
    
    # Calculate the log-likelihood for the single Gaussian model
    # Log-likelihood is a measure of how well the model fits the data (higher is better).
    log_likelihood_single_gaussian = np.sum(norm.logpdf(data, loc=mean, scale=std))

    print("--- Model 1: Single Gaussian ---")
    print(f"Fitting a single Gaussian distribution to the data.")
    print(f"Calculated Mean: {mean:.2f}, Standard Deviation: {std:.2f}")
    print(f"Total Log-Likelihood for Single Gaussian Model: {log_likelihood_single_gaussian:.2f}\n")

    # 3. Fit a Gaussian Mixture Model (GMM) with 2 components
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    # The gmm.score() method returns the average log-likelihood per sample.
    # We multiply by the number of samples to get the total log-likelihood.
    log_likelihood_gmm = gmm.score(data) * len(data)

    print("--- Model 2: Gaussian Mixture Model (GMM) with 2 Components ---")
    print("Fitting a GMM to the same data.")
    gmm_means = gmm.means_.flatten()
    gmm_weights = gmm.weights_.flatten()
    print(f"Found two components with means: {gmm_means[0]:.2f} and {gmm_means[1]:.2f}")
    print(f"Component weights: {gmm_weights[0]:.2f} and {gmm_weights[1]:.2f}")
    print(f"Total Log-Likelihood for GMM: {log_likelihood_gmm:.2f}\n")

    # 4. Conclusion
    print("--- Conclusion ---")
    print("Comparing the total log-likelihood values:")
    print(f"Single Gaussian: {log_likelihood_single_gaussian:.2f}")
    print(f"GMM:             {log_likelihood_gmm:.2f}")
    if log_likelihood_gmm > log_likelihood_single_gaussian:
        print("\nThe GMM has a significantly higher log-likelihood, indicating a much better fit.")
        print("This is because a mixture of Gaussians can model more complex, multi-modal distributions,")
        print("which is a very common characteristic of real-world data.")
    else:
        print("\nThe single Gaussian model performed better, which is unexpected for this data.")

demonstrate_gmm_fit()