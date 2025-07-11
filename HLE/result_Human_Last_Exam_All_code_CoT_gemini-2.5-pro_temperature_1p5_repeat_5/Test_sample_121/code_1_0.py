import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gaussian_mixture():
    """
    Demonstrates the effectiveness of a Gaussian Mixture Model (GMM)
    for modeling multi-modal data compared to a single Gaussian.
    """
    # 1. Create synthetic multi-modal data
    # This simulates a situation with two underlying sub-populations (related to option C)
    np.random.seed(0)
    # Component 1: 300 samples from N(-4, 1.0^2)
    data1 = np.random.normal(loc=-4, scale=1.0, size=300)
    # Component 2: 500 samples from N(4, 1.5^2)
    data2 = np.random.normal(loc=4, scale=1.5, size=500)
    
    # Combine into one dataset
    X = np.concatenate((data1, data2)).reshape(-1, 1)

    # 2. Fit a single Gaussian model
    mu_single = np.mean(X)
    std_single = np.std(X)

    # 3. Fit a Gaussian Mixture Model (GMM) with K=2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(X)

    # Print model parameters
    print("--- Model Fitting Results ---")
    print("\nSingle Gaussian Model:")
    print(f"  - Mean: {mu_single:.4f}")
    print(f"  - Std Dev: {std_single:.4f}")

    print("\nGaussian Mixture Model (GMM) Components:")
    # Sort components by mean for consistent output
    sorted_indices = np.argsort(gmm.means_.flatten())
    weights = gmm.weights_[sorted_indices]
    means = gmm.means_.flatten()[sorted_indices]
    covariances = gmm.covariances_.flatten()[sorted_indices]

    for i in range(gmm.n_components):
        print(f"\nComponent {i+1}:")
        print(f"  - Weight (Mixture Proportion): {weights[i]:.4f}")
        print(f"  - Mean: {means[i]:.4f}")
        # GMM learns variance, so we take sqrt for std dev
        print(f"  - Std Dev (sqrt of covariance): {np.sqrt(covariances[i]):.4f}")
    
    # --- Visualization ---
    plt.figure(figsize=(12, 7))
    
    # Plot histogram of the data
    plt.hist(X, bins=30, density=True, alpha=0.6, color='gray', label='Data Histogram')
    
    # Create a range of x-values for plotting PDFs
    x_plot = np.linspace(X.min(), X.max(), 1000).reshape(-1, 1)

    # Plot the single Gaussian PDF
    pdf_single = norm.pdf(x_plot, mu_single, std_single)
    plt.plot(x_plot, pdf_single, 'r-', linewidth=2, label='Fitted Single Gaussian')
    
    # Plot the GMM PDF
    # The GMM PDF is calculated from its `score_samples` method, which gives log probabilities
    log_pdf_gmm = gmm.score_samples(x_plot)
    plt.plot(x_plot, np.exp(log_pdf_gmm), 'g--', linewidth=3, label='Fitted Gaussian Mixture (GMM)')
    
    plt.title('Comparison of Single Gaussian vs. GMM on Multi-Modal Data')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, linestyle=':')
    
    print("\n--- Conclusion ---")
    print("The plot is being displayed.")
    print("As the visualization shows, the single Gaussian (red line) fails to capture the two-peaked nature of the data.")
    print("The Gaussian Mixture Model (green dashed line) provides an excellent fit, successfully modeling the complex, multi-modal distribution.")
    print("This demonstrates that using a mixture of Gaussians is a powerful technique for real-world data, as stated in option A.")
    
    plt.show()

if __name__ == '__main__':
    demonstrate_gaussian_mixture()