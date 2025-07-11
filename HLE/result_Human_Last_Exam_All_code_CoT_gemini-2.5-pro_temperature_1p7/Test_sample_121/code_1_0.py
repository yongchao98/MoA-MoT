import numpy as np
from sklearn.mixture import GaussianMixture

def solve():
    """
    Illustrates why a Gaussian Mixture Model (GMM) is a better fit for complex data
    than a single Gaussian, supporting option A.
    """
    # For reproducibility
    np.random.seed(0)

    # 1. Generate synthetic bimodal data, simulating a complex real-world distribution.
    # A single Gaussian would be a poor fit for this data.
    data_1 = np.random.normal(loc=-4, scale=1.0, size=300)
    data_2 = np.random.normal(loc=4, scale=1.5, size=700)
    data = np.concatenate([data_1, data_2]).reshape(-1, 1)

    print("--- Model Comparison on Bimodal Data ---")

    # 2. Fit a single Gaussian model (equivalent to a GMM with n_components=1)
    gmm_single = GaussianMixture(n_components=1, random_state=0)
    gmm_single.fit(data)
    bic_single = gmm_single.bic(data)
    print(f"Single Gaussian Model BIC: {bic_single:.2f}")

    # 3. Fit a 2-component GMM
    gmm_mixture = GaussianMixture(n_components=2, random_state=0)
    gmm_mixture.fit(data)
    bic_mixture = gmm_mixture.bic(data)
    print(f"2-Component GMM BIC: {bic_mixture:.2f}")

    # 4. Compare models and print conclusion
    # Lower BIC (Bayesian Information Criterion) indicates a better model fit.
    print("\nConclusion:")
    if bic_mixture < bic_single:
        print("The 2-Component GMM has a significantly lower BIC, indicating a superior fit.")
        print("This demonstrates that using a mixture of Gaussians is a good idea as it can model complex, multi-modal distributions that a single Gaussian cannot capture well.")
    else:
        print("The models have comparable BIC scores in this run.")

    print("\n--- Final Fitted GMM Equation ---")
    print("A GMM's probability density is p(x) = sum(w_i * N(x | mu_i, sigma_i^2))")
    
    # Sort components by mean for consistent output
    sorted_indices = np.argsort(gmm_mixture.means_.flatten())
    weights = gmm_mixture.weights_[sorted_indices]
    means = gmm_mixture.means_.flatten()[sorted_indices]
    covariances = gmm_mixture.covariances_.flatten()[sorted_indices]
    
    # 5. Print the "equation" with the final numbers
    for i, (w, m, c) in enumerate(zip(weights, means, covariances)):
        print(f"\nComponent {i+1}:")
        print(f"  Weight (w{i+1}): {w:.4f}")
        print(f"  Mean (mu{i+1}): {m:.4f}")
        print(f"  Variance (sigma{i+1}^2): {c:.4f}")
        print(f"Equation term {i+1}: {w:.4f} * N(x | {m:.4f}, {c:.4f})")

solve()