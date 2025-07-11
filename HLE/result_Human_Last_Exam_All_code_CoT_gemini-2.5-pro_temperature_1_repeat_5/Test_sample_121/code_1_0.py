import numpy as np
from sklearn.mixture import GaussianMixture

def demonstrate_gmm_fitting():
    """
    This function demonstrates how a Gaussian Mixture Model (GMM) can identify
    the underlying components of a multi-modal distribution.
    """
    # Set a seed for reproducibility
    np.random.seed(0)

    # 1. Create a bimodal dataset by combining two different Gaussian distributions
    # This simulates real-world data that is not a simple, single Gaussian.
    data_1 = np.random.normal(loc=2, scale=1, size=300)  # First component
    data_2 = np.random.normal(loc=9, scale=2, size=700)  # Second component
    data = np.concatenate([data_1, data_2]).reshape(-1, 1)

    # 2. Construct and fit a Gaussian Mixture Model with K=2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # 3. Print the learned parameters for the final GMM equation.
    # The equation for the GMM's probability density is a weighted sum of Gaussian densities:
    # p(x) = w_1 * N(x|μ_1, σ_1^2) + w_2 * N(x|μ_2, σ_2^2)
    # We will now print the learned parameters (w, μ, σ^2) for each component.

    print("A Gaussian Mixture Model was fitted to the bimodal data.")
    print("The parameters for the final probability density equation are:\n")

    for i in range(gmm.n_components):
        weight = gmm.weights_[i]
        mean = gmm.means_[i][0]
        # covariance is variance in 1D
        variance = gmm.covariances_[i][0][0]

        print(f"--- Component {i+1} ---")
        print(f"Weight (w_{i+1}): {weight:.4f}")
        print(f"Mean (μ_{i+1}): {mean:.4f}")
        print(f"Variance (σ^2_{i+1}): {variance:.4f}\n")

    print("As you can see, the GMM successfully identified the two underlying distributions")
    print("(around means of ~2.0 and ~9.0), demonstrating its ability to model complex distributions.")

if __name__ == '__main__':
    demonstrate_gmm_fitting()