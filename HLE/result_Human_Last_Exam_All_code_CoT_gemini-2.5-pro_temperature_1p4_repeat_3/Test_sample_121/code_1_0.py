import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def solve():
    """
    Demonstrates why a Gaussian Mixture Model (GMM) can be better than a single Gaussian.
    """
    # Set a seed for reproducibility
    np.random.seed(0)

    # 1. Create synthetic bimodal data by combining two different Gaussians
    # This represents a complex distribution that a single Gaussian cannot model well.
    data1 = np.random.normal(loc=-4, scale=1.0, size=300)
    data2 = np.random.normal(loc=4, scale=1.5, size=300)
    data = np.concatenate([data1, data2]).reshape(-1, 1)

    print("--- Data Generation ---")
    print("Generated a bimodal dataset by combining samples from two Gaussians:")
    print("Component 1: Mean=-4.0, StdDev=1.0")
    print("Component 2: Mean=4.0, StdDev=1.5\n")

    # 2. Fit a single Gaussian model to the data
    # The "equation" for this model is defined by its mean and standard deviation.
    mu_single = np.mean(data)
    std_single = np.std(data)
    
    print("--- Single Gaussian Model Fit ---")
    print("This model tries to explain all data with one distribution.")
    print(f"Fitted Mean: {mu_single:.4f}")
    print(f"Fitted Standard Deviation: {std_single:.4f}\n")

    # 3. Fit a Gaussian Mixture Model (GMM) with 2 components
    # The "equation" for this model is defined by the weights, means, and variances of its components.
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # Sort components by mean for consistent output
    sorted_indices = np.argsort(gmm.means_.flatten())
    weights_gmm = gmm.weights_[sorted_indices]
    means_gmm = gmm.means_.flatten()[sorted_indices]
    std_devs_gmm = np.sqrt(gmm.covariances_.flatten())[sorted_indices]

    print("--- Gaussian Mixture Model (GMM) Fit ---")
    print("This model explains the data with a mixture of two Gaussians.")
    for i in range(len(weights_gmm)):
        print(f"Component {i+1}:")
        print(f"  - Weight: {weights_gmm[i]:.4f}")
        print(f"  - Mean: {means_gmm[i]:.4f}")
        print(f"  - Standard Deviation: {std_devs_gmm[i]:.4f}")
    print("")

    # 4. Compare the models using log-likelihood (higher is better)
    # The GMM's `score` method returns the average log-likelihood per sample.
    log_likelihood_gmm = gmm.score(data) * len(data)
    # Calculate log-likelihood for the single Gaussian model
    log_likelihood_single = np.sum(norm.logpdf(data, loc=mu_single, scale=std_single))

    print("--- Model Comparison ---")
    print("Log-likelihood measures how well each model fits the data (higher is better).")
    print(f"Total Log-Likelihood (Single Gaussian): {log_likelihood_single:.2f}")
    print(f"Total Log-Likelihood (GMM): {log_likelihood_gmm:.2f}\n")
    
    print("--- Conclusion ---")
    print("The GMM's parameters closely match the original two distributions.")
    print("The GMM's log-likelihood is significantly higher, proving it's a much better fit.")
    print("This demonstrates that a mixture of Gaussians is excellent for modeling complex, multi-modal distributions where a single Gaussian is inaccurate.")


solve()
<<<A>>>