import numpy as np
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt

def demonstrate_gmm_fit():
    """
    This function demonstrates why a Gaussian Mixture Model (GMM) can be
    superior to a single Gaussian for modeling complex data.
    """
    # 1. Generate multi-modal data
    # This data comes from two different "sources" or "clusters"
    np.random.seed(0)
    # Component 1: 300 points centered at 2 with std dev 1.5
    data1 = np.random.normal(loc=2, scale=1.5, size=300)
    # Component 2: 200 points centered at 9 with std dev 1.0
    data2 = np.random.normal(loc=9, scale=1.0, size=200)
    # Combine into a single dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    # 2. Fit a single Gaussian model
    # This is equivalent to a GMM with n_components=1
    single_gaussian = GaussianMixture(n_components=1, random_state=0)
    single_gaussian.fit(data)
    log_likelihood_single = single_gaussian.score(data) # This is the average log-likelihood

    # 3. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)
    log_likelihood_gmm = gmm.score(data) # Average log-likelihood

    # 4. Print results to show the GMM is a better fit
    print("Demonstration: Fitting a Single Gaussian vs. a GMM to multi-modal data.")
    print("-" * 70)
    print(f"Data was generated from two distinct normal distributions.")
    print(f"Mean of the single Gaussian fit: {single_gaussian.means_[0][0]:.2f}")
    print(f"Average log-likelihood for single Gaussian: {log_likelihood_single:.4f}")
    print("-" * 70)
    print(f"Means of the two GMM components: {gmm.means_[0][0]:.2f} and {gmm.means_[1][0]:.2f}")
    print(f"Weights of the two GMM components: {gmm.weights_[0]:.2f} and {gmm.weights_[1]:.2f}")
    print(f"Average log-likelihood for GMM: {log_likelihood_gmm:.4f}")
    print("-" * 70)
    print("Conclusion: The GMM has a significantly higher log-likelihood, indicating a much better fit.")
    print("This illustrates option A: a mixture of Gaussians can model complex, multi-modal distributions effectively.")
    
    # Although not part of the final output, plotting makes it clear.
    # We will generate a plot and save it to a file.
    x_plot = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)
    
    # Get PDF values for the single Gaussian
    log_prob_single = single_gaussian.score_samples(x_plot)
    pdf_single = np.exp(log_prob_single)

    # Get PDF values for the GMM
    log_prob_gmm = gmm.score_samples(x_plot)
    pdf_gmm = np.exp(log_prob_gmm)

    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=30, density=True, alpha=0.5, label='Data Histogram')
    plt.plot(x_plot, pdf_single, color='red', linestyle='--', linewidth=2, label='Single Gaussian Fit')
    plt.plot(x_plot, pdf_gmm, color='green', linewidth=2, label='Gaussian Mixture (K=2) Fit')
    plt.title('Single Gaussian vs. Gaussian Mixture Model')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True)
    plt.savefig('gmm_comparison.png')
    print("\nA plot illustrating this comparison has been saved to 'gmm_comparison.png'")

demonstrate_gmm_fit()