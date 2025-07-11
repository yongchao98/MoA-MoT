import numpy as np
from sklearn.mixture import GaussianMixture

def evaluate_emission_density_models():
    """
    Demonstrates that a Gaussian Mixture Model (GMM) provides a better fit 
    for complex, multi-modal data compared to a single Gaussian.
    """
    # Seed for reproducibility
    np.random.seed(0)

    # 1. Create synthetic data from two different Gaussian distributions
    # This creates a bimodal (two-peaked) dataset
    data1 = np.random.normal(loc=-2, scale=0.8, size=300)
    data2 = np.random.normal(loc=3, scale=1.2, size=300)
    data = np.concatenate([data1, data2]).reshape(-1, 1)

    print("Evaluating models on a bimodal dataset...")
    print("-" * 40)

    # 2. Fit a model with a single Gaussian component (K=1)
    gmm_single = GaussianMixture(n_components=1, random_state=0)
    gmm_single.fit(data)
    bic_single = gmm_single.bic(data)
    
    # 3. Fit a model with two Gaussian components (K=2)
    gmm_mixture = GaussianMixture(n_components=2, random_state=0)
    gmm_mixture.fit(data)
    bic_mixture = gmm_mixture.bic(data)

    # 4. Print the results for comparison. A lower BIC is better.
    print(f"Model with 1 Gaussian component:")
    print(f"  - BIC Score: {bic_single:.2f}")

    print(f"\nModel with 2 Gaussian components (a mixture):")
    print(f"  - BIC Score: {bic_mixture:.2f}")
    
    print("\n--- Model Parameters for the Mixture (K=2) ---")
    weights = gmm_mixture.weights_
    means = gmm_mixture.means_
    covariances = gmm_mixture.covariances_

    print("The final probability density is a sum of weighted Gaussians.")
    for i in range(len(weights)):
        print(f"\nComponent {i+1}:")
        print(f"  - The weight is: {weights[i]:.4f}")
        print(f"  - The mean (mu) is: {means[i][0]:.4f}")
        print(f"  - The variance (sigma^2) is: {covariances[i][0][0]:.4f}")

    print("-" * 40)
    print("\nConclusion: The model with a mixture of Gaussians has a much lower BIC score,")
    print("indicating it's a significantly better fit for the multi-modal data.")


evaluate_emission_density_models()