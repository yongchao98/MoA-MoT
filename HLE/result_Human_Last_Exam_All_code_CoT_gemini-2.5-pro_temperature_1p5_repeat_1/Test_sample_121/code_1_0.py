import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Demonstrates that a Gaussian Mixture Model (GMM) provides a better fit
    than a single Gaussian for non-Gaussian (bimodal) data.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Generate synthetic bimodal data
    # This data comes from two different Gaussian distributions, so a single
    # Gaussian will not be able to model it accurately.
    data_component1 = np.random.normal(-4, 1.2, 300)  # 300 points from the first mode
    data_component2 = np.random.normal(3, 1.5, 700)   # 700 points from the second mode
    data = np.concatenate([data_component1, data_component2]).reshape(-1, 1)

    print("--- A Demonstration of Model Fitting for Bimodal Data ---")
    print(f"Generated a dataset of {len(data)} points from two distinct Gaussians.\n")

    # 2. Fit a single Gaussian distribution to the data
    # We use norm.fit() from scipy to find the best-fit mu and sigma
    mu, std = norm.fit(data)

    # Calculate the log-likelihood for the single Gaussian model
    log_likelihood_sg = np.sum(norm.logpdf(data, loc=mu, scale=std))
    # Calculate AIC. k=2 parameters (mu, std). AIC = 2*k - 2*log-likelihood
    aic_sg = 2 * 2 - 2 * log_likelihood_sg

    print("--- Model 1: Single Gaussian Fit ---")
    print("The 'equation' for this model's probability density is defined by:")
    print(f"  Mean (μ): {mu:.4f}")
    print(f"  Std Dev (σ): {std:.4f}")
    print(f"\nModel Quality Score (lower is better):")
    print(f"  Akaike Information Criterion (AIC): {aic_sg:.2f}\n")


    # 3. Fit a Gaussian Mixture Model (with 2 components) to the data
    # We know the ground truth has 2 components, so we set n_components=2
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    # The GMM object can calculate its own AIC
    aic_gmm = gmm.aic(data)

    print("--- Model 2: Gaussian Mixture Model (GMM) Fit ---")
    print("The 'equation' for this model's probability density is a weighted sum of components:")
    # Loop through the components and print their parameters
    weights = gmm.weights_
    means = gmm.means_
    covariances = gmm.covariances_
    for i in range(len(weights)):
        print(f"  Component {i+1}:")
        print(f"    Weight: {weights[i]:.4f}")
        print(f"    Mean (μ): {means[i][0]:.4f}")
        print(f"    Std Dev (σ): {np.sqrt(covariances[i][0][0]):.4f}")

    print(f"\nModel Quality Score (lower is better):")
    print(f"  Akaike Information Criterion (AIC): {aic_gmm:.2f}\n")

    # 4. Conclusion
    print("--- Conclusion ---")
    if aic_gmm < aic_sg:
        print("The GMM has a significantly lower AIC, indicating it is a much better model for this data.")
        print("This supports the idea that using a mixture of Gaussians is effective for complex, non-Gaussian data.")
    else:
        print("In this case, the single Gaussian performed better, which is unexpected.")

if __name__ == '__main__':
    demonstrate_gmm_fit()