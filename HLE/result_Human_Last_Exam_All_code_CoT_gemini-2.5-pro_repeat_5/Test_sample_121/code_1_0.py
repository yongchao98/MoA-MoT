import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    This function demonstrates the superiority of a Gaussian Mixture Model (GMM)
    for fitting complex, multi-modal data compared to a single Gaussian distribution.
    """
    # Step 1: Generate synthetic data from two modes to simulate a complex distribution.
    # This represents real-world data which often doesn't fit a simple bell curve.
    np.random.seed(42)
    data_mode1 = np.random.normal(loc=-4, scale=1.5, size=300)
    data_mode2 = np.random.normal(loc=5, scale=2.0, size=700)
    data = np.concatenate([data_mode1, data_mode2]).reshape(-1, 1)

    # Step 2: Fit a single Gaussian distribution.
    # This model assumes the data has only one central peak.
    single_gaussian_mean = data.mean()
    single_gaussian_std = data.std()

    # Step 3: Fit a Gaussian Mixture Model with K=2 components.
    # This model can identify and represent the two underlying clusters in the data.
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    # Step 4: Prepare for plotting by creating a range of x values.
    x_range = np.linspace(-12, 12, 1000).reshape(-1, 1)

    # Calculate the PDF for the single Gaussian model.
    pdf_single = norm.pdf(x_range, single_gaussian_mean, single_gaussian_std)

    # Calculate the PDF for the GMM. `score_samples` gives log-probabilities.
    log_pdf_gmm = gmm.score_samples(x_range)
    pdf_gmm = np.exp(log_pdf_gmm)

    # --- Output and Explanation ---
    print("--- Analysis of Model Fits ---")
    print("This script demonstrates why using a mixture of Gaussians is a good idea.")
    print("\nObservation from the plot (if displayed):")
    print(" - The data histogram clearly shows two distinct peaks (it's bimodal).")
    print(" - The 'Single Gaussian Fit' (red line) tries to cover all the data with one curve, failing to capture the bimodal nature.")
    print(" - The 'Gaussian Mixture Model Fit' (green line) accurately models both peaks, providing a much better representation of the data's true distribution.")
    
    print("\nConclusion:")
    print("This directly supports the idea that a mixture of Gaussians is effective because it can model more complex, multi-modal distributions often found in real-world data.")

    # --- Plotting ---
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=50, density=True, alpha=0.7, color='skyblue', label='Data Histogram')
    plt.plot(x_range, pdf_single, 'r-', linewidth=2, label='Single Gaussian Fit')
    plt.plot(x_range, pdf_gmm, 'g--', linewidth=3, label='Gaussian Mixture Model Fit (K=2)')
    plt.title('Single Gaussian vs. Gaussian Mixture Model on Multi-modal Data')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    demonstrate_gmm_fit()