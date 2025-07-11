import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    This function demonstrates why a Gaussian Mixture Model (GMM) is a good
    choice for modeling complex, multi-modal data compared to a single Gaussian.
    """
    # 1. Generate synthetic 'real-world' data from two different groups
    # This creates a bimodal (two-peaked) distribution.
    np.random.seed(42)
    # Data from the first subpopulation (e.g., one operational mode)
    data1 = np.random.normal(-4, 1.2, 300)
    # Data from the second subpopulation (e.g., another operational mode)
    data2 = np.random.normal(5, 1.8, 700)
    # Combine into a single dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    print("--- Demonstrating Model Fits ---")
    print("We have created a dataset with two distinct groups to simulate complex real-world data.")
    print("Let's see how well a single Gaussian and a Gaussian Mixture Model (GMM) can fit it.")

    # 2. Fit a single Gaussian distribution
    # This model assumes the data has only one central peak.
    mu_single, std_single = norm.fit(data)
    print("\n--- Single Gaussian Fit ---")
    print(f"This model incorrectly assumes all data comes from one source.")
    print(f"Fitted Mean: {mu_single:.4f}, Fitted Standard Deviation: {std_single:.4f}")
    print("A single Gaussian model is a poor fit because it averages out the two distinct peaks, failing to capture the underlying structure.")

    # 3. Fit a Gaussian Mixture Model (GMM) with K=2 components
    # This model assumes the data is a mix of K (here, K=2) different Gaussian distributions.
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    print("\n--- Gaussian Mixture Model (GMM) Fit ---")
    print("This model correctly assumes the data can come from multiple sources (components).")
    print("The final model is a weighted sum of the components, representing the final equation for the density.")
    
    # 4. Output the parameters (the numbers in the final equation)
    # The GMM density equation is: p(x) = weight_1 * Gaussian_1(x) + weight_2 * Gaussian_2(x)
    print("\nFinal GMM Equation Parameters:")
    for i in range(gmm.n_components):
        weight = gmm.weights_[i]
        mean = gmm.means_[i][0]
        # GMM gives covariance, so we take its square root for standard deviation
        std_dev = np.sqrt(gmm.covariances_[i][0][0])
        
        print(f"\n--- Component {i+1} ---")
        print(f"This component represents one of the underlying groups in the data.")
        print(f"Weight (how much of the data belongs to this group): {weight:.4f}")
        print(f"Mean (center of this group): {mean:.4f}")
        print(f"Standard Deviation (spread of this group): {std_dev:.4f}")

    print("\n--- Conclusion ---")
    print("The GMM successfully identified the two underlying distributions in our complex data.")
    print("This demonstrates that using a mixture of Gaussians is a powerful technique for modeling complex, multi-modal distributions often found in real-world data, where a single Gaussian would be inaccurate.")
    print("This directly validates option A.")

if __name__ == '__main__':
    demonstrate_gmm_fit()