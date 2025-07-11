import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Demonstrates that a Gaussian Mixture Model (GMM) can fit a bimodal
    distribution better than a single Gaussian.
    """
    # 1. Generate synthetic bimodal data
    np.random.seed(0)
    # Data from first Gaussian component
    data1 = np.random.normal(loc=-4, scale=1.5, size=300)
    # Data from second Gaussian component
    data2 = np.random.normal(loc=5, scale=2, size=700)
    # Combine into one dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    # 2. Fit a single Gaussian to the data
    mean_single = np.mean(data)
    std_single = np.std(data)

    # 3. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # Extract GMM parameters
    weights = gmm.weights_
    means = gmm.means_
    covariances = gmm.covariances_
    stds = np.sqrt(covariances)

    # 4. Print the parameters for the final GMM equation
    print("--- Single Gaussian Fit ---")
    print(f"This model assumes a single equation: N(x | mean, std^2)")
    print(f"Mean: {mean_single:.4f}")
    print(f"Standard Deviation: {std_single:.4f}\n")

    print("--- Gaussian Mixture Model (GMM) Fit ---")
    print("This model uses a mixture equation: p(x) = sum(weight_k * N(x | mean_k, std_k^2))")
    print(f"Found {gmm.n_components} components:")
    for i in range(gmm.n_components):
        print(f"\n--- Component {i+1} ---")
        print(f"Equation Term: {weights[i]:.4f} * N(x | {means[i][0]:.4f}, {stds[i][0]**2:.4f})")
        print(f"Weight: {weights[i]:.4f}")
        print(f"Mean: {means[i][0]:.4f}")
        print(f"Standard Deviation: {stds[i][0]:.4f}")

    print("\nAs you can see, the GMM correctly identifies the two underlying distributions")
    print("from the data, whereas the single Gaussian finds an inaccurate average.")

if __name__ == '__main__':
    demonstrate_gmm_fit()
