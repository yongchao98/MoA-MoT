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

    print("--- Step 1: Generating a Bimodal (Two-Peaked) Dataset ---")
    # Generate data from two different Gaussian distributions
    data1 = np.random.normal(loc=-4, scale=1.5, size=300)
    data2 = np.random.normal(loc=5, scale=2.0, size=700)
    
    # Combine them into a single bimodal dataset
    data = np.concatenate([data1, data2]).reshape(-1, 1)
    print(f"Generated a dataset with {data.shape[0]} points.")
    print("The data is intentionally non-Gaussian, with two distinct peaks.\n")

    print("--- Step 2: Fitting a Single Gaussian Model ---")
    # Calculate parameters for a single Gaussian fit
    mean_sg = np.mean(data)
    std_sg = np.std(data)
    
    # Calculate the log-likelihood for the single Gaussian model
    # Log-likelihood is the sum of the log probabilities of each data point
    log_likelihood_sg = np.sum(norm.logpdf(data, loc=mean_sg, scale=std_sg))
    
    print(f"Single Gaussian Model Parameters:")
    print(f"  - Mean = {mean_sg:.4f}")
    print(f"  - Std. Dev. = {std_sg:.4f}")
    print(f"Total Log-Likelihood for Single Gaussian Model: {log_likelihood_sg:.4f}\n")

    print("--- Step 3: Fitting a Gaussian Mixture Model (GMM) ---")
    # Fit a GMM with 2 components
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)
    
    # The GMM object provides a method to score the data (mean log-likelihood)
    # We multiply by the number of data points to get the total log-likelihood
    log_likelihood_gmm = gmm.score(data) * data.shape[0]
    
    print("Gaussian Mixture Model (2 Components) Parameters:")
    for i in range(gmm.n_components):
        print(f"  - Component {i+1}:")
        print(f"    - Weight = {gmm.weights_[i]:.4f}")
        print(f"    - Mean   = {gmm.means_[i][0]:.4f}")
        print(f"    - Std. Dev. (from covariance) = {np.sqrt(gmm.covariances_[i][0][0]):.4f}")

    print(f"\nTotal Log-Likelihood for GMM: {log_likelihood_gmm:.4f}\n")

    print("--- Step 4: Conclusion ---")
    print("Comparing the log-likelihoods:")
    print(f"  Single Gaussian: {log_likelihood_sg:.4f}")
    print(f"  GMM:             {log_likelihood_gmm:.4f}")
    
    if log_likelihood_gmm > log_likelihood_sg:
        print("\nThe GMM has a significantly higher log-likelihood, indicating a much better fit.")
        print("This demonstrates that a mixture of Gaussians can model complex, multi-modal distributions")
        print("more accurately than a single Gaussian, supporting option A.")
    else:
        print("\nThe single Gaussian model had a better or equal fit, which is unexpected for this data.")

if __name__ == '__main__':
    demonstrate_gmm_fit()

<<<A>>>