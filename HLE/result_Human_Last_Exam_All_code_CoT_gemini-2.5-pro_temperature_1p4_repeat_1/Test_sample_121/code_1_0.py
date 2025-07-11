import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Generates bimodal data and fits both a single Gaussian and a GMM
    to demonstrate the GMM's ability to model complex distributions.
    """
    # 1. Generate synthetic bimodal data
    np.random.seed(42)
    # Component 1
    data1 = np.random.normal(loc=-4, scale=1, size=300)
    # Component 2
    data2 = np.random.normal(loc=4, scale=1.5, size=700)
    # Combine into one dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    # 2. Fit a single Gaussian distribution
    mu_single, std_single = norm.fit(data)
    
    # 3. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)
    
    # --- Output the parameters (the "equation") ---
    print("--- Single Gaussian Fit ---")
    print(f"This model assumes all data comes from one distribution.")
    # The equation for a single Gaussian is defined by its mean and standard deviation.
    print(f"Mean (μ): {mu_single:.4f}")
    print(f"Standard Deviation (σ): {std_single:.4f}")
    print("\n")

    print("--- Gaussian Mixture Model (GMM) Fit ---")
    print(f"This model assumes the data is a mix of {gmm.n_components} distributions.")
    # The equation for a GMM is a weighted sum of Gaussians.
    for i in range(gmm.n_components):
        print(f"Component {i+1}:")
        print(f"  Weight: {gmm.weights_[i]:.4f}")
        print(f"  Mean (μ): {gmm.means_[i][0]:.4f}")
        # Covariance is variance for 1D data. Std dev is sqrt of variance.
        print(f"  Standard Deviation (σ): {np.sqrt(gmm.covariances_[i][0][0]):.4f}")
    
    # 4. Plot the results for visualization
    plt.figure(figsize=(12, 6))
    
    # Plot the histogram of the data
    plt.hist(data, bins=50, density=True, alpha=0.6, label='Data Histogram')
    
    # Plot the single Gaussian PDF
    x = np.linspace(data.min(), data.max(), 1000)
    pdf_single = norm.pdf(x, mu_single, std_single)
    plt.plot(x, pdf_single, 'r-', linewidth=2, label='Single Gaussian Fit (Poor Fit)')
    
    # Plot the GMM PDF
    logprob = gmm.score_samples(x.reshape(-1, 1))
    pdf_gmm = np.exp(logprob)
    plt.plot(x, pdf_gmm, 'g-', linewidth=2, label='Gaussian Mixture Fit (Good Fit)')

    plt.title('Comparison of Single Gaussian vs. GMM Fit on Bimodal Data')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

# Run the demonstration
demonstrate_gmm_fit()