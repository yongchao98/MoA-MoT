import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Generates multi-modal data and compares the fit of a single Gaussian
    versus a Gaussian Mixture Model (GMM).
    """
    # 1. Generate synthetic multi-modal data
    # This simulates real-world data that isn't from a single Gaussian distribution.
    np.random.seed(0)
    # Data from the first mode (cluster)
    data1 = np.random.normal(loc=-4, scale=1.0, size=300)
    # Data from the second mode (cluster)
    data2 = np.random.normal(loc=5, scale=1.5, size=500)
    # Data from a third, wider mode (could represent outliers or another sub-population)
    data3 = np.random.normal(loc=0, scale=3, size=100)
    
    # Combine into one dataset
    data = np.concatenate([data1, data2, data3]).reshape(-1, 1)

    # 2. Fit a single Gaussian distribution
    mu, std = norm.fit(data)

    # 3. Fit a Gaussian Mixture Model with 3 components
    gmm = GaussianMixture(n_components=3, random_state=0)
    gmm.fit(data)

    # 4. Prepare data for plotting the PDFs
    x_plot = np.linspace(-10, 10, 1000).reshape(-1, 1)
    
    # PDF for the single Gaussian
    pdf_single_gaussian = norm.pdf(x_plot, mu, std)
    
    # PDF for the GMM (weighted sum of component Gaussians)
    log_prob = gmm.score_samples(x_plot)
    pdf_gmm = np.exp(log_prob)

    # 5. Plot the results
    plt.figure(figsize=(12, 7))
    plt.hist(data, bins=50, density=True, alpha=0.6, color='gray', label='Observed Data Histogram')
    
    plt.plot(x_plot, pdf_single_gaussian, 'r-', linewidth=2, label='Single Gaussian Fit')
    plt.plot(x_plot, pdf_gmm, 'b--', linewidth=2, label='Gaussian Mixture Model (GMM) Fit')
    
    plt.title('Comparison of Single Gaussian vs. GMM Fit')
    plt.xlabel('Data Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Print a conclusion based on the plot
    print("Demonstration Results:")
    print("-----------------------")
    print(f"Single Gaussian Fit: Mean={mu:.2f}, StdDev={std:.2f}")
    print("The single Gaussian (red line) fails to capture the multi-modal nature of the data.")
    print("\nGaussian Mixture Model Fit (3 components):")
    for i in range(gmm.n_components):
        print(f"  Component {i+1}: Weight={gmm.weights_[i]:.2f}, Mean={gmm.means_[i][0]:.2f}, StdDev={np.sqrt(gmm.covariances_[i][0][0]):.2f}")
    print("The GMM (blue dashed line) provides a much better fit by modeling each peak separately.")
    print("\nConclusion: Using a mixture of Gaussians is a good idea because it can model complex distributions, unlike a single Gaussian. This supports option A.")
    
    # The plot will be displayed if running in an interactive environment.
    # To save it, you could uncomment the following line:
    # plt.savefig("gmm_fit_comparison.png")
    plt.show()

# Run the demonstration
demonstrate_gmm_fit()