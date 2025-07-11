import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def solve_and_explain():
    """
    This function demonstrates why a Gaussian Mixture Model (GMM) can be
    superior to a single Gaussian distribution for modeling real-world data.
    """
    # Step 1: Create a bimodal dataset by combining two different normal distributions.
    # A single Gaussian distribution will be a poor fit for this kind of data.
    np.random.seed(0)
    data1 = np.random.normal(loc=-4, scale=1.5, size=300)
    data2 = np.random.normal(loc=5, scale=2.0, size=700)
    bimodal_data = np.concatenate([data1, data2]).reshape(-1, 1)

    print("--- Data Analysis ---")
    print(f"Generated a dataset with {len(bimodal_data)} points, which is intentionally bimodal.")
    print("\n")

    # Step 2: Fit a single Gaussian distribution to the data.
    # We get its mean and standard deviation.
    mu, std = norm.fit(bimodal_data)
    print("--- Single Gaussian Fit ---")
    print("A single Gaussian distribution attempts to model this data with one peak.")
    print(f"Fitted Mean (μ): {mu:.4f}")
    print(f"Fitted Standard Deviation (σ): {std:.4f}")
    print("\n")

    # Step 3: Fit a Gaussian Mixture Model (GMM) with K=2 components.
    # This is more appropriate for bimodal data.
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(bimodal_data)

    # Step 4: Print the parameters of the fitted GMM.
    # This corresponds to the components of the "final equation". The emission
    # density P(observation) is a weighted sum of these Gaussian components.
    # P(x) = weight_1 * N(x | mean_1, cov_1) + weight_2 * N(x | mean_2, cov_2)
    print("--- Gaussian Mixture Model (GMM) Fit ---")
    print("A GMM can model the two underlying clusters in the data.")
    print("The final 'equation' for the emission density is a sum of the following Gaussians:\n")

    for i in range(gmm.n_components):
        weight = gmm.weights_[i]
        mean = gmm.means_[i][0]
        covariance = gmm.covariances_[i][0][0]
        
        print(f"Component {i+1}:")
        print(f"  Weight (π_{i+1}): {weight:.4f}")
        print(f"  Mean (μ_{i+1}):   {mean:.4f}")
        print(f"  Variance (σ²_{i+1}): {covariance:.4f}")
        print("-" * 20)

    print("\nConclusion: The GMM successfully identifies two distinct components close to")
    print("our original distributions (means around -4 and 5), providing a far more")
    print("accurate model of the data's complex, multi-modal distribution.")


# Run the demonstration
solve_and_explain()