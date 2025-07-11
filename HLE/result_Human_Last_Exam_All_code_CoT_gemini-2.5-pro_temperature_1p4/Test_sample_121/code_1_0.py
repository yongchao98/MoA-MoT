import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import matplotlib.pyplot as plt

def solve():
    """
    This function demonstrates why a Gaussian Mixture Model (GMM) can be a better
    fit for complex data than a single Gaussian distribution.
    """
    # Set a random seed for reproducibility
    np.random.seed(0)

    # 1. Generate synthetic bimodal data
    # Data from the first mode
    data1 = np.random.normal(loc=-4, scale=1.5, size=300)
    # Data from the second mode
    data2 = np.random.normal(loc=5, scale=2.0, size=700)
    # Combine into one dataset
    data = np.concatenate((data1, data2)).reshape(-1, 1)

    # 2. Fit a single Gaussian distribution
    single_mean, single_std = norm.fit(data)

    # 3. Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=0)
    gmm.fit(data)

    # Extract GMM parameters
    weights = gmm.weights_
    means = gmm.means_
    covariances = gmm.covariances_
    stds = np.sqrt(covariances)

    # 4. Output the parameters and equations for the fitted models
    print("--- Fitting a Single Gaussian Distribution ---")
    print(f"This model assumes all data comes from one source.")
    print(f"The equation for the probability density function (PDF) is: N(x | μ, σ²)")
    print("Fitted Parameters:")
    print(f"  μ (mean)     = {single_mean:.4f}")
    print(f"  σ (std dev)  = {single_std:.4f}")
    print("\nFinal Equation:")
    print(f"PDF(x) = N(x | {single_mean:.4f}, {single_std**2:.4f})")
    print("-" * 45)

    print("\n--- Fitting a Gaussian Mixture Model (GMM) ---")
    print("This model assumes the data is a weighted sum of multiple Gaussian sources.")
    print("The equation for the PDF is: Σ [weight_k * N(x | μ_k, σ²_k)]")
    print("Fitted Parameters:")
    for i in range(gmm.n_components):
        print(f"\nComponent {i+1}:")
        print(f"  Weight (π_k) = {weights[i]:.4f}")
        print(f"  μ_k (mean)   = {means[i][0]:.4f}")
        print(f"  σ_k (std dev)= {stds[i][0]:.4f}")
    
    print("\nFinal Equation:")
    eq_parts = []
    for i in range(gmm.n_components):
        part = f"{weights[i]:.4f} * N(x | {means[i][0]:.4f}, {covariances[i][0][0]:.4f})"
        eq_parts.append(part)
    print(f"PDF(x) = {' + '.join(eq_parts)}")
    print("-" * 45)


    # 5. Visualize the results
    x_plot = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

    # PDF for the single Gaussian
    single_gauss_pdf = norm.pdf(x_plot, single_mean, single_std)

    # PDF for the GMM
    # The score_samples method gives the log-probability density, so we exponentiate it
    log_gmm_pdf = gmm.score_samples(x_plot)
    gmm_pdf = np.exp(log_gmm_pdf)
    
    plt.figure(figsize=(12, 7))
    plt.hist(data, bins=30, density=True, alpha=0.6, color='skyblue', label='Data Histogram')
    plt.plot(x_plot, single_gauss_pdf, color='red', linewidth=2, linestyle='--', label='Single Gaussian Fit')
    plt.plot(x_plot, gmm_pdf, color='green', linewidth=2, label='Gaussian Mixture Model Fit (K=2)')
    
    plt.title('GMM vs. Single Gaussian Fit on Bimodal Data')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    print("\nA plot has been generated to visually compare the fits.")
    print("Close the plot window to end the script.")
    plt.show()

# Run the solution
solve()