import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Generates bimodal data and fits both a single Gaussian and a GMM
    to demonstrate the GMM's superior ability to model complex distributions.
    """
    # 1. Generate synthetic data from a mixture of two Gaussians
    np.random.seed(0)
    # Component 1
    data1 = np.random.normal(loc=-4, scale=1, size=300)
    # Component 2
    data2 = np.random.normal(loc=4, scale=1.5, size=700)
    # Combine into one dataset
    data = np.concatenate([data1, data2]).reshape(-1, 1)

    # 2. Fit a single Gaussian distribution
    mu_single = np.mean(data)
    sigma_single = np.std(data)

    # 3. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # 4. Prepare for plotting
    x_plot = np.linspace(-10, 10, 1000).reshape(-1, 1)

    # PDF for the single Gaussian
    pdf_single = norm.pdf(x_plot, mu_single, sigma_single)

    # PDF for the GMM
    # The GMM pdf is a weighted sum of the component pdfs
    pdf_gmm = np.exp(gmm.score_samples(x_plot))

    # Print model parameters
    print("--- Single Gaussian Fit ---")
    print(f"Mean: {mu_single:.2f}")
    print(f"Standard Deviation: {sigma_single:.2f}")
    print("\n--- Gaussian Mixture Model (GMM) Fit ---")
    print(f"Component Weights: {gmm.weights_[0]:.2f}, {gmm.weights_[1]:.2f}")
    print(f"Component Means: {gmm.means_[0][0]:.2f}, {gmm.means_[1][0]:.2f}")
    print(f"Component Covariances (Variances): {gmm.covariances_[0][0][0]:.2f}, {gmm.covariances_[1][0][0]:.2f}")


    # 5. Plot the results
    plt.figure(figsize=(12, 6))
    plt.hist(data, bins=50, density=True, alpha=0.6, color='gray', label='Data Histogram')
    plt.plot(x_plot, pdf_single, 'r-', linewidth=2, label='Single Gaussian Fit')
    plt.plot(x_plot, pdf_gmm, 'g--', linewidth=2, label='Gaussian Mixture Model Fit (K=2)')
    plt.title('Fitting a Bimodal Distribution')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

# Run the demonstration
demonstrate_gmm_fit()