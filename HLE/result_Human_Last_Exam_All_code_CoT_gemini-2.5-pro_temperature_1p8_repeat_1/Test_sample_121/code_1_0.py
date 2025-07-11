import numpy as np
from sklearn.mixture import GaussianMixture

def solve_task():
    """
    Demonstrates the effectiveness of a Gaussian Mixture Model (GMM)
    for multi-modal data compared to a single Gaussian model.
    """
    # --- Step 1: Generate a bimodal dataset ---
    # This simulates real-world data that isn't a simple bell curve.
    np.random.seed(42)
    # Create data from two distinct distributions
    data_component1 = np.random.normal(loc=-5, scale=1.5, size=300)
    data_component2 = np.random.normal(loc=5, scale=2, size=700)
    data = np.concatenate((data_component1, data_component2)).reshape(-1, 1)

    print("--- Analysis of Fitting Multi-Modal Data ---")
    print("The data was generated from two distinct Gaussian distributions.")
    print(f"Ground Truth Component 1: Mean=-5.0, Weight=30%")
    print(f"Ground Truth Component 2: Mean= 5.0, Weight=70%\n")


    # --- Step 2: Fit a single Gaussian model (equivalent to a GMM with one component) ---
    single_gaussian_model = GaussianMixture(n_components=1, random_state=42)
    single_gaussian_model.fit(data)

    print("--- Result 1: Fitting a Single Gaussian Model ---")
    print("The final equation for the probability density is: N(x | mean, variance)")
    sg_mean = single_gaussian_model.means_[0][0]
    sg_cov = single_gaussian_model.covariances_[0][0]
    print(f"Fitted Mean     = {sg_mean:.4f}")
    print(f"Fitted Variance = {sg_cov:.4f}")
    print("Observation: The single Gaussian finds a mean somewhere between the two actual data clusters, failing to accurately represent either.\n")


    # --- Step 3: Fit a Gaussian Mixture Model with two components ---
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    print("--- Result 2: Fitting a Gaussian Mixture Model (GMM) ---")
    print("The final equation is: weight_1 * N(x | mean_1, var_1) + weight_2 * N(x | mean_2, var_2)")
    # Sort the components by their means for consistent ordering
    indices = np.argsort(gmm.means_.flatten())
    weights = gmm.weights_[indices]
    means = gmm.means_.flatten()[indices]
    covariances = gmm.covariances_.flatten()[indices]

    print("Fitted Parameters for each component in the equation:")
    for i in range(gmm.n_components):
        print(f"  Component {i+1}:")
        print(f"    Weight   (weight_{i+1}) = {weights[i]:.4f}")
        print(f"    Mean     (mean_{i+1})   = {means[i]:.4f}")
        print(f"    Variance (var_{i+1})    = {covariances[i]:.4f}")

    print("\nObservation: The GMM successfully identifies the two underlying clusters in the data, with its fitted means and weights closely matching the ground truth. This illustrates its ability to model complex, multi-modal distributions.")

solve_task()
<<<A>>>