import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import lasso_path
from sklearn.datasets import make_regression

def demonstrate_lasso_equivalence_nuance():
    """
    This function illustrates the relationship between the two forms of Lasso
    regression and highlights the nuances that make the equivalence "not always
    strictly true".
    """
    # 1. Generate a synthetic dataset.
    # We choose p > n (n_features > n_samples) to create a situation
    # where the solution path can have interesting properties.
    n_samples, n_features = 50, 80
    X, y = make_regression(n_samples=n_samples, n_features=n_features,
                           n_informative=10, noise=5.0, random_state=42)

    # 2. Compute the full Lasso path
    # We use a fine grid of many alpha (lambda) values to see the path clearly.
    # The path is computed from a large lambda down to a small one.
    alphas, coefs, _ = lasso_path(X, y, n_alphas=1000, eps=1e-6)

    # 3. Calculate t(lambda), the L1 norm of coefficients for each lambda
    t_values = np.sum(np.abs(coefs), axis=0)

    # 4. Illustrate the relationship with a specific example
    print("--- Illustrating the Lasso Equivalence ---")
    # Pick a point on the path
    idx_example = 500
    lambda_example = alphas[idx_example]
    beta_example = coefs[:, idx_example]
    t_example = t_values[idx_example]

    print("The two Lasso formulations are:")
    print("1. Constrained: argmin(RSS) s.t. sum(|beta_j|) <= t")
    print("2. Penalized:   argmin(RSS + lambda * sum(|beta_j|))")
    print("\nThe equivalence means a solution for one is a solution for the other.")
    print("\nFor example:")
    print(f"For a penalty lambda = {lambda_example:.4f}, we find a solution where the L1 norm is t = {t_example:.4f}.")
    # Printing the first 5 coefficients of the solution vector for brevity
    print(f"The first 5 coefficients of the solution vector are: {np.round(beta_example[:5], 2)}")
    print("This means this specific solution vector solves both the penalized problem with the given lambda,")
    print("and the constrained problem with the resulting t as the budget.")


    # 5. Check for nuances that make the equivalence "not always strictly true"
    # The equivalence is "not strict" if the mapping between lambda and t is not
    # a clean one-to-one function. This happens if t(lambda) is flat.
    # We check for this by looking at the differences between consecutive t_values.
    # Note: We reverse the arrays because lasso_path returns alphas in decreasing order.
    t_values_sorted = t_values[::-1]
    alphas_sorted = alphas[::-1]
    
    # Using np.isclose to handle floating point inaccuracies
    diff_t = np.diff(t_values_sorted)
    flat_spots_indices = np.where(np.isclose(diff_t, 0))[0]

    print("\n--- Checking for 'Strict' Equivalence Nuances ---")
    if len(flat_spots_indices) > 0:
        flat_spot_t = t_values_sorted[flat_spots_indices[0]]
        flat_spot_lambda_start = alphas_sorted[flat_spots_indices[0]]
        flat_spot_lambda_end = alphas_sorted[flat_spots_indices[0] + 1]

        print(f"We found 'flat spots' in the solution path.")
        print(f"For instance, the L1 norm t remains constant at ~{flat_spot_t:.4f} while lambda increases")
        print(f"from {flat_spot_lambda_start:.4f} to {flat_spot_lambda_end:.4f}.")
        print("\nThis means for a single 't', there can be a range of corresponding 'lambda' values.")
        print("This breaks the one-to-one mapping and is why the equivalence is often described as")
        print("'true in general, but not always strictly true'.")
    else:
        print("\nIn this simulation, the L1 norm 't' was a strictly monotonic function of lambda.")
        print("This represents the 'well-behaved' general case where the equivalence appears strict.")

    # 6. Plotting the results for visualization
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Coefficient paths
    ax1.plot(np.log10(alphas), coefs.T)
    ax1.set_xlabel('log10(lambda)')
    ax1.set_ylabel('Coefficients (beta_j)')
    ax1.set_title('Lasso Coefficient Paths')
    ax1.set_xlim(np.log10(alphas.min()), np.log10(alphas.max()))

    # Plot 2: t vs. lambda
    ax2.plot(alphas, t_values, lw=2)
    ax2.set_xlabel('Penalty (lambda)')
    ax2.set_ylabel('L1 Norm (t = sum(|beta_j|))')
    ax2.set_title('L1 Norm vs. Penalty')
    ax2.set_xscale('log')
    ax2.set_xlim(alphas.min(), alphas.max())
    if len(flat_spots_indices) > 0:
        ax2.axhline(y=flat_spot_t, color='r', linestyle='--', lw=1.5, label=f'Flat spot at t={flat_spot_t:.2f}')
        ax2.legend()
    
    fig.suptitle("Visualizing the Lasso Path and t-lambda Relationship", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# Run the demonstration
demonstrate_lasso_equivalence_nuance()
