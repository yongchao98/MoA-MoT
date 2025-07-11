import numpy as np
from sklearn.mixture import GaussianMixture

def demonstrate_gmm_fit():
    """
    Demonstrates that a Gaussian Mixture Model (GMM) provides a better fit
    for multi-modal data compared to a single Gaussian distribution.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Generate a multi-modal dataset
    # Create data from two different Gaussian distributions and combine them
    data_1 = np.random.normal(loc=-4, scale=1.0, size=300)
    data_2 = np.random.normal(loc=4, scale=1.5, size=300)
    data = np.concatenate([data_1, data_2]).reshape(-1, 1)

    # 2. Fit a single Gaussian model (GMM with 1 component)
    gmm_single = GaussianMixture(n_components=1, random_state=42)
    gmm_single.fit(data)
    bic_single = gmm_single.bic(data)

    # 3. Fit a Gaussian Mixture Model (GMM with 2 components)
    gmm_mixture = GaussianMixture(n_components=2, random_state=42)
    gmm_mixture.fit(data)
    bic_mixture = gmm_mixture.bic(data)

    # 4. Print and explain the results
    print("Demonstration: Fitting multi-modal data with different models.")
    print("-" * 60)
    print("We created a dataset with two distinct peaks (a bimodal distribution).")
    print("We will now compare how well a single Gaussian vs. a mixture of two Gaussians can model this data.")
    print("We use the Bayesian Information Criterion (BIC) for comparison - a lower score is better.")
    print("\n--- Model Comparison ---")

    print(f"\nModel 1: Single Gaussian (n_components=1)")
    print(f"BIC Score: {bic_single:.2f}")

    print(f"\nModel 2: Gaussian Mixture (n_components=2)")
    print(f"BIC Score: {bic_mixture:.2f}")

    print("\n--- Conclusion ---")
    if bic_mixture < bic_single:
        print("The Gaussian Mixture Model has a significantly lower BIC score.")
        print("This indicates it is a much better fit for the complex, multi-modal data.")
        print("This illustrates the principle that using a mixture of Gaussians is a good idea to model complex distributions often found in real-world data.")
    else:
        print("The single Gaussian model performed better, which is unexpected for this bimodal data.")

if __name__ == '__main__':
    demonstrate_gmm_fit()