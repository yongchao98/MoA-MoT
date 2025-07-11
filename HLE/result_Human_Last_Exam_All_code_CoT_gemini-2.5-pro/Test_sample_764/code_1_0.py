import numpy as np
from sklearn.linear_model import Lasso
from scipy.optimize import minimize

def solve_lasso_equivalence_puzzle():
    """
    This function demonstrates that the equivalence between the constrained and
    penalized forms of Lasso is not always strictly true, using a p > n counterexample.
    """
    # Step 1 & 2: Set up a p > n scenario where a perfect fit (RSS=0) is possible.
    # n=2 observations, p=3 predictors.
    X = np.array([[1., 2., 0.],
                  [2., 1., 1.]])
    y = np.array([3., 5.])

    # The general solution for X*beta = y can be parameterized by a single variable 'c':
    # beta = (3 - 2c, c, 3c - 1)
    
    # Step 3: Find the solution beta* that minimizes the L1 norm among all perfect-fit solutions.
    # This is the solution the penalized form will find as lambda -> 0.
    def l1_norm_of_solution(c):
        """Calculates the L1 norm for a solution parameterized by c."""
        c = c[0]
        beta = np.array([3 - 2*c, c, 3*c - 1])
        return np.linalg.norm(beta, 1)

    # Find the value of 'c' that minimizes the L1 norm. Analytically, it's c=1/3.
    res = minimize(l1_norm_of_solution, x0=[0.0])
    c_star = res.x[0]
    beta_star = np.array([3 - 2*c_star, c_star, 3*c_star - 1])
    t_star = l1_norm_of_solution([c_star])

    print("--- The Minimum L1-Norm Solution (beta*) ---")
    print(f"The minimum L1-norm solution (beta*) is: ({beta_star[0]:.4f}, {beta_star[1]:.4f}, {beta_star[2]:.4f})")
    print(f"The minimum L1 norm (t*) is: {t_star:.4f}")
    print(f"RSS for beta*: {np.sum((y - X @ beta_star)**2):.4f}\n")

    # Step 4: Define another perfect-fit solution beta' with a larger L1 norm.
    # This beta' is a valid solution for the constrained problem with budget t' > t*.
    # Let's pick c = 0.5, which gives a L1 norm of 3.
    c_prime = 0.5
    beta_prime = np.array([3 - 2*c_prime, c_prime, 3*c_prime - 1])
    t_prime = np.linalg.norm(beta_prime, 1)

    print("--- Another Perfect-Fit Solution (beta') ---")
    print(f"The alternative solution (beta') is: ({beta_prime[0]:.4f}, {beta_prime[1]:.4f}, {beta_prime[2]:.4f})")
    print(f"The L1 norm of beta' (t') is: {t_prime:.4f}")
    print(f"RSS for beta': {np.sum((y - X @ beta_prime)**2):.4f}\n")

    # Step 5 & 6: Show that the penalized form (sklearn's Lasso) finds beta*, not beta'.
    print("--- Analysis of Equivalence ---")
    print(f"The solution beta' is a valid solution to the constrained problem for a budget t' = {t_prime:.4f}.")
    print("However, for the penalized problem, any lambda > 0 will select the solution with the minimum L1 norm.")
    
    # Use a very small lambda (alpha) to approximate the lambda -> 0+ case.
    lasso_solver = Lasso(alpha=1e-8, fit_intercept=False, tol=1e-8, max_iter=10000)
    lasso_solver.fit(X, y)
    beta_lasso = lasso_solver.coef_

    print("\n--- Penalized Lasso Solver Result ---")
    print(f"For a very small lambda (alpha={lasso_solver.alpha}), the computed coefficients are: ({beta_lasso[0]:.4f}, {beta_lasso[1]:.4f}, {beta_lasso[2]:.4f})")
    print("This is effectively beta*, not beta'.")
    print("\nConclusion: A solution from the constrained problem (beta') cannot be found by the penalized")
    print("problem for any lambda > 0. Thus, the equivalence is not strictly true.")
    
    # Final output with numbers as requested
    print(f"\nThe final equation showing the L1 norm of the penalized Lasso solution ({np.linalg.norm(beta_lasso, 1):.4f}) is approximately t* ({t_star:.4f}) and not t' ({t_prime:.4f}), which illustrates the breakdown in strict equivalence.")

if __name__ == '__main__':
    solve_lasso_equivalence_puzzle()