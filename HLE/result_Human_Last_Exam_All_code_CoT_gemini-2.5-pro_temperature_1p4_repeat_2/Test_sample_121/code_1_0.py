import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    This script demonstrates why using a mixture of Gaussians is a good idea for
    modeling complex data, as stated in option A. We will create a bimodal
    dataset that a single Gaussian cannot model accurately, but a Gaussian Mixture
    Model (GMM) can.
    """
    # 1. Generate a complex, bimodal dataset
    np.random.seed(42)
    # Component 1: 300 samples from N(-4, 1.0^2)
    data1 = np.random.normal(loc=-4, scale=1.0, size=300)
    # Component 2: 500 samples from N(4, 1.5^2)
    data2 = np.random.normal(loc=4, scale=1.5, size=500)
    # Combine into a single dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    # 2. Fit a single Gaussian distribution
    # The "equation" for a single Gaussian is defined by its mean and standard deviation
    single_gaussian_mean = data.mean()
    single_gaussian_std = data.std()
    
    # 3. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)
    
    # The "equation" for a GMM is defined by the weights, means, and covariances
    # of each component.
    gmm_weights = gmm.weights_
    gmm_means = gmm.means_
    gmm_covariances = gmm.covariances_

    # 4. Print the "equations" (parameters) for each model
    print("--- Model Fitting Results ---")
    print("\n[Single Gaussian Model]")
    print(f"This model assumes all data comes from one distribution.")
    print("Equation Parameters:")
    print(f"  mean                = {single_gaussian_mean:.4f}")
    print(f"  standard_deviation  = {single_gaussian_std:.4f}")

    print("\n[Gaussian Mixture Model (GMM) with K=2]")
    print(f"This model correctly assumes the data is a mix of two distributions.")
    print("Equation Parameters:")
    for i in range(len(gmm_weights)):
        print(f"  Component {i+1}:")
        print(f"    Weight      = {gmm_weights[i]:.4f}")
        print(f"    Mean        = {gmm_means[i][0]:.4f}")
        print(f"    Std. Dev.   = {np.sqrt(gmm_covariances[i][0][0]):.4f}")
        
    # 5. Plot the results for visualization
    plt.figure(figsize=(12, 7))
    x_axis = np.linspace(-10, 10, 1000).reshape(-1, 1)

    # Plot histogram of the actual data
    plt.hist(data, bins=30, density=True, alpha=0.6, color='skyblue', label='Actual Data Distribution')
    
    # Plot the PDF of the single Gaussian fit
    single_gaussian_pdf = norm.pdf(x_axis, single_gaussian_mean, single_gaussian_std)
    plt.plot(x_axis, single_gaussian_pdf, 'r--', linewidth=2, label='Single Gaussian Fit (Poor)')
    
    # Plot the PDF of the GMM fit
    logprob = gmm.score_samples(x_axis)
    gmm_pdf = np.exp(logprob)
    plt.plot(x_axis, gmm_pdf, 'g-', linewidth=2, label='Gaussian Mixture Fit (Good)')

    plt.title('GMM vs. Single Gaussian on Multi-Modal Data')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    
    print("\n--- Conclusion ---")
    print("The plot will be displayed. It visually confirms that the GMM provides a much better fit to the complex,")
    print("multi-modal data, supporting the reasoning in Option A.")

    plt.show()

if __name__ == '__main__':
    demonstrate_gmm_fit()