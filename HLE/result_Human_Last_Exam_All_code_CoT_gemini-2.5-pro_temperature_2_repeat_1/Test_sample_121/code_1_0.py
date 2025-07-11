import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def demonstrate_gmm_fit():
    """
    This function demonstrates that a Gaussian Mixture Model (GMM) provides a
    better fit than a single Gaussian for complex, multi-modal data.
    """
    # --- Plan ---
    # 1. Generate bimodal data that a single Gaussian cannot model well.
    # 2. Fit a single Gaussian model and calculate its total log-likelihood.
    # 3. Fit a Gaussian Mixture Model (GMM) and calculate its total log-likelihood.
    # 4. Compare the log-likelihoods to numerically show that the GMM is a better fit,
    #    thus supporting the reasoning in option A.

    # Set a seed for reproducibility of the random data
    np.random.seed(42)

    # 1. Generate synthetic bimodal data
    # This data has two distinct peaks, making it non-Gaussian overall.
    print("Step 1: Generating a bimodal dataset with 1000 observations.")
    data1 = np.random.normal(loc=-4, scale=1, size=500)
    data2 = np.random.normal(loc=4, scale=1.5, size=500)
    data = np.concatenate((data1, data2)).reshape(-1, 1)
    print("Dataset created.\n")

    # 2. Fit a single Gaussian model
    print("--- Fitting a Single Gaussian Model ---")
    # A single Gaussian is defined by its mean and standard deviation.
    # The 'equation' for its log-likelihood is: sum(log(N(x_i | mu, sigma^2)))
    mean_single = np.mean(data)
    std_single = np.std(data)

    # Calculate the total log-likelihood for the single Gaussian model.
    # A higher log-likelihood indicates a better model fit.
    log_likelihood_single = np.sum(norm.logpdf(data, loc=mean_single, scale=std_single))

    print(f"Fitted a single Gaussian with Mean = {mean_single:.4f} and Std Dev = {std_single:.4f}")
    print("The final 'equation' value for the log-likelihood of this model is:")
    print(f"Total Log-Likelihood: {log_likelihood_single:.4f}\n")


    # 3. Fit a Gaussian Mixture Model (GMM) with K=2 components
    print("--- Fitting a Gaussian Mixture Model (GMM) ---")
    # The 'equation' for GMM log-likelihood is: sum(log(sum_k(w_k * N(x_i | mu_k, sigma_k^2))))
    # We set n_components=2 because we know our synthetic data has two modes.
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(data)

    # The 'score' method gives the average log-likelihood per sample.
    # We multiply by the number of samples to get the total log-likelihood.
    log_likelihood_gmm = gmm.score(data) * len(data)

    # We can inspect the parameters the GMM found for each component
    weights = gmm.weights_
    means = gmm.means_
    stds = np.sqrt(gmm.covariances_).flatten()
    
    print(f"Fitted a GMM with 2 components:")
    print(f"  Component 1: Weight={weights[0]:.2f}, Mean={means[0][0]:.4f}, Std Dev={stds[0]:.4f}")
    print(f"  Component 2: Weight={weights[1]:.2f}, Mean={means[1][0]:.4f}, Std Dev={stds[1]:.4f}")
    print("\nThe final 'equation' value for the log-likelihood of this model is:")
    print(f"Total Log-Likelihood: {log_likelihood_gmm:.4f}\n")

    # 4. Compare the results and conclude
    print("--- Comparison ---")
    print("To choose the best model, we compare their log-likelihood values.")
    print(f"Single Gaussian Log-Likelihood: {log_likelihood_single:.4f}")
    print(f"GMM Log-Likelihood:             {log_likelihood_gmm:.4f}")
    
    if log_likelihood_gmm > log_likelihood_single:
        improvement = log_likelihood_gmm - log_likelihood_single
        print(f"\nThe GMM provides a substantially better fit (improvement of {improvement:.2f}).")
        print("This numerically demonstrates that using a mixture of Gaussians can model complex distributions more accurately, validating option A.")
    else:
        print("\nUnexpectedly, the single Gaussian provided a better or equal fit.")

demonstrate_gmm_fit()