import numpy as np
from sklearn.linear_model import Lasso

def demonstrate_lasso_equivalence():
    """
    Demonstrates that the equivalence between the penalized and constrained
    LASSO formulations is not always strictly one-to-one.
    """
    # Set up a dataset (X, y) where y is orthogonal to the columns of X.
    # This ensures the Ordinary Least Squares (OLS) solution is beta=0.
    np.random.seed(42)
    n_samples, n_features = 20, 5
    X = np.random.randn(n_samples, n_features)
    
    # Create a y vector that is orthogonal to the column space of X
    # by subtracting its projection onto that space.
    y_random = np.random.randn(n_samples)
    y_proj_X = X @ np.linalg.pinv(X) @ y_random
    y = y_random - y_proj_X

    # For this data, the LASSO solution should be beta=0 for any lambda >= 0.
    # Let's verify this for different values of lambda (alpha in scikit-learn).
    
    lambda_values = [0.1, 1.0, 10.0]
    
    print("Investigating the relationship between lambda and t in LASSO formulations.\n")

    for lam in lambda_values:
        # Solve the penalized LASSO problem: argmin RSS + lambda * sum(|beta|)
        # The 'alpha' parameter in sklearn's Lasso is our lambda.
        # We don't fit an intercept since y is centered w.r.t X's projection.
        lasso_model = Lasso(alpha=lam, fit_intercept=False)
        lasso_model.fit(X, y)
        
        # Get the solution beta_hat
        beta_hat = lasso_model.coef_
        
        # Calculate the corresponding t value: t = sum(|beta_hat|)
        t_value = np.sum(np.abs(beta_hat))
        
        print(f"For lambda = {lam}:")
        print(f"  Coefficients (beta): {beta_hat}")
        # Final Equation Output: We show the values for lambda, beta_j, and t
        print(f"  Resulting sum of absolute coefficients (t) = " + " + ".join([f"|{b:.4f}|" for b in beta_hat]) + f" = {t_value:.4f}")
        print("-" * 20)

    print("\nConclusion:")
    print("As demonstrated, multiple values of the penalty parameter lambda (0.1, 1.0, 10.0)")
    print("all lead to the same solution where all beta coefficients are zero.")
    print("This results in the same value for the constraint parameter t = sum(|beta_j|) = 0.")
    print("\nThis shows a many-to-one mapping from lambda to t, meaning the equivalence")
    print("between the two LASSO formulations is not a strict (one-to-one) correspondence.")
    print("The statement is true in that an equivalence exists, but it is not 'strictly' true.")

demonstrate_lasso_equivalence()
<<<K>>>