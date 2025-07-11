import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Generates bimodal data and demonstrates the superiority of a
    Gaussian Mixture Model (GMM) fit over a single Gaussian fit.
    """
    # --- Step 1: Generate multi-modal data ---
    # This data represents observations from a single state that has a complex,
    # non-Gaussian distribution. We create it by combining two different Gaussian distributions.
    np.random.seed(0)
    # Data from the first mode/cluster
    data1 = np.random.normal(-5, 1.5, 300)
    # Data from the second mode/cluster
    data2 = np.random.normal(5, 2.0, 700)
    # Combine them into a single dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    # --- Step 2: Fit a single Gaussian (the inaccurate model) ---
    # A single Gaussian will try to find the overall mean and standard deviation,
    # failing to capture the two separate peaks.
    mu_single, std_single = norm.fit(data)

    # --- Step 3: Fit a Gaussian Mixture Model (the better model) ---
    # A GMM with 2 components can model the two underlying distributions separately.
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # --- Step 4: Print the parameters of the fitted models ("final equations") ---
    print("--- Model Fitting Results ---")
    print("\n[Inaccurate Model] Single Gaussian Parameters:")
    print(f"Equation: N(x | μ={mu_single:.2f}, σ^2={std_single**2:.2f})")
    print(f"μ (mean) = {mu_single:.2f}")
    print(f"σ (std dev) = {std_single:.2f}")


    print("\n[Accurate Model] Gaussian Mixture Model Parameters:")
    print("Equation: Σ [ w_i * N(x | μ_i, σ_i^2) ] for i=1 to 2")
    # Loop through each number in the final equation and print it
    for i in range(gmm.n_components):
        weight = gmm.weights_[i]
        mean = gmm.means_[i][0]
        variance = gmm.covariances_[i][0][0]
        print(f"  Component {i+1}:")
        print(f"    w_{i+1} (weight)      = {weight:.2f}")
        print(f"    μ_{i+1} (mean)        = {mean:.2f}")
        print(f"    σ^2_{i+1} (variance)    = {variance:.2f}")

    # --- Step 5: Visualize the results for comparison ---
    x_plot = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)
    pdf_single_gaussian = norm.pdf(x_plot, mu_single, std_single)
    pdf_gmm = np.exp(gmm.score_samples(x_plot))

    plt.figure(figsize=(12, 6))
    plt.hist(data, bins=50, density=True, alpha=0.6, label='Real-world Data Histogram')
    plt.plot(x_plot, pdf_single_gaussian, 'r-', linewidth=2, label='Single Gaussian Fit (Inaccurate)')
    plt.plot(x_plot, pdf_gmm, 'g--', linewidth=3, label='Gaussian Mixture Fit (Accurate)')
    plt.title('Fitting a Multi-Modal Distribution')
    plt.xlabel('Observation Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    demonstrate_gmm_fit()