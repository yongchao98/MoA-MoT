import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.mixture import GaussianMixture

def solve_and_plot():
    """
    This function demonstrates why a Gaussian Mixture Model (GMM) is superior
    to a single Gaussian for modeling multi-modal data.
    """
    # 1. Generate synthetic bimodal data
    # Create two separate distributions and combine them
    np.random.seed(42) # for reproducibility
    data1 = np.random.normal(loc=-4, scale=1.5, size=300)
    data2 = np.random.normal(loc=5, scale=2, size=700)
    data = np.concatenate([data1, data2]).reshape(-1, 1)

    # 2. Fit a single Gaussian distribution
    mu, std = norm.fit(data)

    # 3. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    # 4. Prepare data for plotting the PDFs
    x_plot = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

    # Get PDF for the single Gaussian
    pdf_single_gaussian = norm.pdf(x_plot, mu, std)

    # Get PDF for the GMM
    # The GMM's score_samples gives log probabilities, so we exponentiate
    log_pdf_gmm = gmm.score_samples(x_plot)
    pdf_gmm = np.exp(log_pdf_gmm)

    # 5. Plot the results
    plt.figure(figsize=(12, 7))
    plt.hist(data, bins=50, density=True, alpha=0.6, color='skyblue', label='Real Data Distribution')
    plt.plot(x_plot, pdf_single_gaussian, 'r-', linewidth=2, label='Single Gaussian Fit')
    plt.plot(x_plot, pdf_gmm, 'k--', linewidth=3, label='Gaussian Mixture Model Fit (K=2)')

    plt.title('Fitting a Multi-Modal Distribution', fontsize=16)
    plt.xlabel('Value', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    
    print("Plot generated successfully. Observe how the GMM (black dashed line) accurately captures")
    print("the two peaks in the data, while the single Gaussian (red solid line) fails to do so,")
    print("placing its single peak in the valley between the real data clusters.")
    print("\nThis visually confirms that a mixture of Gaussians can model more complex,")
    print("multi-modal distributions where a single Gaussian is inaccurate.")

    plt.show()

# Run the demonstration
solve_and_plot()