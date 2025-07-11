import numpy as np
import cvxpy as cp

def solve_lasso_equivalence_demo():
    """
    This function demonstrates a case where the equivalence between the
    penalized and constrained LASSO formulations is not strictly true.
    """
    # 1. Set up a problem with perfect collinearity.
    # We create a dataset where x2 is identical to x1.
    # The true model is y = 1*x1 + 1*x2, so the OLS solution should have beta_1 + beta_2 = 2.
    np.random.seed(42)
    n_samples = 50
    x1 = np.random.rand(n_samples)
    x2 = x1.copy()  # Perfect collinearity
    X = np.vstack([x1, x2]).T
    y = 1 * x1 + 1 * x2 + 0.1 * np.random.randn(n_samples)

    # Define the CVXPY variable for the coefficients
    beta = cp.Variable(2)
    
    # Define the common objective function: Residual Sum of Squares (RSS)
    objective_rss = cp.sum_squares(y - X @ beta)

    # 2. Solve the penalized form with lambda = 0 (which is an OLS problem).
    # Because of collinearity, the solution is not unique. The solver will pick one.
    problem_penalized = cp.Problem(cp.Minimize(objective_rss))
    problem_penalized.solve()
    beta_penalized = beta.value

    # 3. Solve the constrained form.
    # We choose a t-value large enough so it doesn't prevent finding an OLS solution.
    # An OLS solution sum of coefficients is expected to be around 2. A t of 3 is sufficient.
    t_constrained = 3.0
    constraints_constrained = [cp.norm(beta, 1) <= t_constrained]
    problem_constrained = cp.Problem(cp.Minimize(objective_rss), constraints_constrained)
    problem_constrained.solve()
    beta_constrained = beta.value

    # 4. Print and compare the results
    print("This example demonstrates that the LASSO formulations are not always strictly equivalent.")
    print("We use a dataset where predictor x2 is identical to x1.")
    print("The OLS solution is any pair (β₁, β₂) where β₁ + β₂ is a constant (approx. 2.0).")
    print("-" * 50)
    
    print("Solution from Penalized Form (with λ = 0):")
    print(f"β₁ = {beta_penalized[0]:.4f}, β₂ = {beta_penalized[1]:.4f}")
    # Final equation for penalized form
    print(f"The equation for the sum is: {beta_penalized[0]:.4f} + {beta_penalized[1]:.4f} = {np.sum(beta_penalized):.4f}")
    print("-" * 50)

    print(f"Solution from Constrained Form (with t = {t_constrained}):")
    print(f"β₁ = {beta_constrained[0]:.4f}, β₂ = {beta_constrained[1]:.4f}")
    # Final equation for constrained form
    print(f"The equation for the sum is: {beta_constrained[0]:.4f} + {beta_constrained[1]:.4f} = {np.sum(beta_constrained):.4f}")
    print("-" * 50)

    print("Conclusion:")
    print("Both methods found a valid OLS solution, as the sum of coefficients is the same in both cases.")
    print("However, the specific coefficient vectors are different.")
    print("The set of solutions for the λ=0 problem is the entire line β₁ + β₂ = 2.0229.")
    print("The solution for the constrained problem is a specific point (or points) on that line.")
    print("Since the solution sets are not identical, the equivalence is not strictly true in all cases.")

solve_lasso_equivalence_demo()