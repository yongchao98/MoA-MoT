import numpy as np
from sklearn.linear_model import Lasso
from sklearn.datasets import make_regression

def solve_lasso_equivalence():
    """
    This function demonstrates the relationship between the constrained and
    penalized forms of Lasso regression to answer the user's question.
    """
    # Step 1: Generate some sample data
    X, y, true_coef = make_regression(n_samples=50, n_features=10, n_informative=5, noise=15, coef=True, random_state=42)
    n_samples = X.shape[0]

    # Step 2: Determine the theoretical lambda_max (alpha_max in sklearn's terms)
    # In sklearn, the objective is (1/(2*n_samples)) * RSS + alpha * L1_norm.
    # The stationarity condition at beta=0 implies that alpha_max is the smallest alpha
    # for which the solution is the zero vector. This occurs when alpha >= max(|X^T*y / n_samples|).
    alpha_max = np.max(np.abs(np.dot(X.T, y))) / n_samples

    print("Demonstrating the non-unique mapping from t to lambda:")
    print("We will show that for a single t value (t=0), multiple lambda (alpha) values give the same result.")
    print("-" * 70)
    print(f"The smallest alpha that should produce all-zero coefficients is alpha_max = {alpha_max:.4f}")
    print("-" * 70)

    # Step 3: Choose multiple alpha values >= alpha_max
    alphas_to_test = [alpha_max, alpha_max + 1.0, alpha_max * 2]

    for i, alpha_val in enumerate(alphas_to_test):
        # Fit the Lasso model (penalized form)
        lasso = Lasso(alpha=alpha_val)
        lasso.fit(X, y)

        # The solution from the penalized form
        beta_hat = lasso.coef_
        
        # The corresponding 't' from the constrained form
        t_hat = np.sum(np.abs(beta_hat))

        print(f"Case {i+1}:")
        print(f"Chosen lambda (alpha) = {alpha_val:.4f}")
        print(f"Resulting coefficients beta_hat = {np.round(beta_hat, 4)}")
        print(f"This corresponds to t = sum(|beta_hat|) = {t_hat:.4f}")
        print("\nThis demonstrates that for the constrained problem with t = 0,")
        print(f"the solution (beta_hat = 0) corresponds to multiple lambda values (e.g., {alpha_val:.4f}).")
        print("-" * 70)

    print("Conclusion:")
    print("The equivalence between the two Lasso formulations is 'True in general' due to convex duality.")
    print("However, it is 'not always strictly true' because the mapping between the parameters t and lambda is not always one-to-one, as shown by the t=0 case.")

solve_lasso_equivalence()