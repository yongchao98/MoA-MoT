import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.mixture import GaussianMixture

def demonstrate_gmm_fit():
    """
    This function demonstrates why a Gaussian Mixture Model (GMM) is superior to a
    single Gaussian for modeling multi-modal data.
    """
    # Step 1: Create a synthetic bimodal dataset
    # This data has two distinct peaks, which a single Gaussian cannot capture.
    np.random.seed(42)
    data_cluster1 = np.random.normal(loc=-3, scale=0.8, size=300)
    data_cluster2 = np.random.normal(loc=4, scale=1.5, size=700)
    data = np.concatenate([data_cluster1, data_cluster2]).reshape(-1, 1)

    # Step 2: Fit a single Gaussian to the data
    # This will find the single best-fit normal distribution.
    mu_single, std_single = norm.fit(data)

    # Step 3: Fit a Gaussian Mixture Model with 2 components
    # We choose K=2 because we know the data has two underlying sources.
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    # --- Output and Visualization ---
    
    # Create a range of x-values for plotting the PDFs
    x_axis = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)
    
    # Plot the data histogram
    plt.hist(data, bins=40, density=True, alpha=0.6, color='skyblue', label='Data Histogram')
    
    # Plot the single Gaussian PDF
    pdf_single_gaussian = norm.pdf(x_axis, mu_single, std_single)
    plt.plot(x_axis, pdf_single_gaussian, 'r-', lw=2, label='Single Gaussian Fit')
    
    # Plot the GMM PDF
    logprob = gmm.score_samples(x_axis)
    pdf_gmm = np.exp(logprob)
    plt.plot(x_axis, pdf_gmm, 'g--', lw=2, label='Gaussian Mixture Fit (K=2)')
    
    plt.title("Single Gaussian vs. Gaussian Mixture Model Fit")
    plt.xlabel("Data Value")
    plt.ylabel("Probability Density")
    plt.legend()
    
    # Save the plot instead of showing it directly
    output_filename = "gmm_vs_gaussian_fit.png"
    plt.savefig(output_filename)

    print(f"A plot comparing the fits has been saved to '{output_filename}'.")
    print("It visually demonstrates that the GMM provides a much more accurate model for the bimodal data.\n")
    
    # Step 4: Print the parameters ("the numbers in the final equation") for the GMM
    # The equation for the GMM's probability density function (PDF) is:
    # p(x) = weight1 * N(x | mean1, std1^2) + weight2 * N(x | mean2, std2^2)
    print("="*60)
    print("Fitted Gaussian Mixture Model (GMM) Parameters")
    print("="*60)
    print("The final equation for the emission density p(x) is composed of these numbers:")
    
    if gmm.converged_:
        # Sort components by their means for consistent ordering
        sorted_indices = np.argsort(gmm.means_.flatten())
        weights = gmm.weights_[sorted_indices]
        means = gmm.means_.flatten()[sorted_indices]
        covariances = gmm.covariances_.flatten()[sorted_indices]

        for i in range(gmm.n_components):
            mean = means[i]
            std_dev = np.sqrt(covariances[i])
            weight = weights[i]
            
            print(f"\n--- Component {i+1} ---")
            print(f"  Weight (w{i+1}): {weight:.4f}")
            print(f"  Mean   (μ{i+1}): {mean:.4f}")
            print(f"  Std Dev (σ{i+1}): {std_dev:.4f}")
    else:
        print("\nThe GMM fitting process did not converge.")

if __name__ == '__main__':
    demonstrate_gmm_fit()