import numpy as np
from sklearn.mixture import GaussianMixture

def demonstrate_gmm_fit():
    """
    Generates bimodal data and fits both a single Gaussian and a GMM
    to demonstrate the GMM's superior ability to model complex distributions.
    """
    # Set a seed for reproducibility
    np.random.seed(0)

    # 1. Generate synthetic bimodal data
    # Create data from two different Gaussian distributions
    n_samples_1 = 500
    n_samples_2 = 500
    
    # Distribution 1: Mean = 2, Std Dev = 1.5
    mean1, cov1 = 2, 1.5**2
    data1 = np.random.normal(mean1, np.sqrt(cov1), n_samples_1)

    # Distribution 2: Mean = 9, Std Dev = 1.0
    mean2, cov2 = 9, 1.0**2
    data2 = np.random.normal(mean2, np.sqrt(cov2), n_samples_2)

    # Combine them into one dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)
    np.random.shuffle(data)

    print("--- Data Generation ---")
    print(f"Generated {n_samples_1} points from N(mean={mean1}, variance={cov1:.2f})")
    print(f"Generated {n_samples_2} points from N(mean={mean2}, variance={cov2:.2f})\n")


    # 2. Fit a single Gaussian distribution
    # The equation is N(mean, variance)
    single_gaussian_mean = np.mean(data)
    single_gaussian_var = np.var(data)
    
    print("--- Single Gaussian Fit ---")
    print("This model assumes all data comes from one distribution.")
    print("Equation: N(mean, variance)")
    print(f"Mean: {single_gaussian_mean:.4f}")
    print(f"Variance: {single_gaussian_var:.4f}\n")


    # 3. Fit a Gaussian Mixture Model with 2 components
    # The equation is w1*N(m1, v1) + w2*N(m2, v2)
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # Sort components by mean for consistent output
    sorted_indices = np.argsort(gmm.means_.flatten())
    weights = gmm.weights_[sorted_indices]
    means = gmm.means_.flatten()[sorted_indices]
    variances = gmm.covariances_.flatten()[sorted_indices]

    print("--- Gaussian Mixture Model (GMM) Fit ---")
    print("This model assumes the data is a mix of 2 distributions.")
    print("Equation: weight1 * N(mean1, variance1) + weight2 * N(mean2, variance2)\n")

    for i in range(len(weights)):
        print(f"Component {i+1}:")
        print(f"  Weight: {weights[i]:.4f}")
        print(f"  Mean: {means[i]:.4f}")
        print(f"  Variance: {variances[i]:.4f}")
        print("-" * 15)

    print("\nConclusion: The GMM successfully identifies the two original distributions,")
    print("while the single Gaussian finds a poor average between them.")

if __name__ == '__main__':
    demonstrate_gmm_fit()