import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    Generates bimodal data and fits both a single Gaussian and a GMM
    to demonstrate the GMM's ability to model complex distributions.
    """
    # 1. Generate a bimodal dataset
    # This data comes from two distinct clusters/subpopulations
    np.random.seed(0)
    data_1 = np.random.normal(loc=-3, scale=1.2, size=300)
    data_2 = np.random.normal(loc=4, scale=1.5, size=700)
    data = np.concatenate([data_1, data_2]).reshape(-1, 1)

    # 2. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # 3. Fit a single Gaussian for comparison
    mean_single = data.mean()
    std_single = data.std()

    # 4. Prepare for plotting
    x_plot = np.linspace(-10, 10, 1000).reshape(-1, 1)

    # Calculate Probability Density Function (PDF) for the single Gaussian
    pdf_single = norm.pdf(x_plot, mean_single, std_single)

    # Calculate PDF for the GMM
    # gmm.score_samples gives log-probabilities, so we exponentiate
    log_prob_gmm = gmm.score_samples(x_plot)
    pdf_gmm = np.exp(log_prob_gmm)

    # 5. Print results and create the plot
    print("Demonstration of GMM vs. Single Gaussian Fit\n")
    print(f"Single Gaussian Fit: Mean={mean_single:.2f}, Std Dev={std_single:.2f}")
    print("\nGaussian Mixture Model Fit (2 Components):")
    for i in range(gmm.n_components):
        print(f"  Component {i+1}: Weight={gmm.weights_[i]:.2f}, Mean={gmm.means_[i][0]:.2f}, Std Dev={np.sqrt(gmm.covariances_[i][0][0]):.2f}")
    print("\nThe plot visually confirms the GMM's superior fit to the complex data.")

    plt.figure(figsize=(12, 7))
    plt.hist(data, bins=60, density=True, alpha=0.6, color='lightblue', label='Data Histogram')
    plt.plot(x_plot, pdf_single, 'r-', linewidth=2, label='Single Gaussian Fit')
    plt.plot(x_plot, pdf_gmm, 'g--', linewidth=2.5, label='Gaussian Mixture Model Fit (K=2)')
    plt.title('GMM vs. Single Gaussian on Bimodal Data')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

if __name__ == '__main__':
    demonstrate_gmm_fit()