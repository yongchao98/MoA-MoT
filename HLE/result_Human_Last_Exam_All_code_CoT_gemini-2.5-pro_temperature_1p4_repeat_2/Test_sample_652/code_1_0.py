import numpy as np
from scipy.optimize import brentq

def solve_and_calculate_elasticity():
    """
    This function solves the economic model to find the elasticity of optimal 
    job search intensity (q) with respect to the probability of unemployment (p).
    """

    # Define the value of p around which to calculate the elasticity.
    p_val = 0.2

    # Step 1: Define the implicit function F(q, p) = 0 which determines the optimal q*.
    # This function is derived from the first-order conditions of the government's 
    # utility maximization problem for the worker.
    # The equation is: (1-q)*exp(2q) - p*(1-q)^2 - (1-p)*(1.5 - q) = 0
    def F(q, p):
        return (1 - q) * np.exp(2 * q) - p * (1 - q)**2 - (1 - p) * (1.5 - q)

    # Step 2: Numerically solve for q* when p = 0.2.
    # We are looking for a root of the function F(q, 0.2) in the interval (0, 1).
    # brentq is an efficient and robust root-finding algorithm.
    try:
        # Search for the root in a slightly smaller interval to avoid issues at the boundaries.
        q_star = brentq(lambda q: F(q, p_val), 0.0001, 0.9999)
    except ValueError:
        print("Error: Could not find the optimal q* in the specified interval.")
        return

    # Step 3: Define the partial derivatives of F needed for implicit differentiation.
    # We need dF/dp and dF/dq to calculate dq/dp.
    def dF_dp(q, p):
        """Calculates the partial derivative of F with respect to p."""
        return -(1 - q)**2 + (1.5 - q)

    def dF_dq(q, p):
        """Calculates the partial derivative of F with respect to q."""
        term1 = np.exp(2 * q) * (1 - 2 * q)
        term2 = 2 * p * (1 - q)
        term3 = 1 - p
        return term1 + term2 + term3

    # Step 4: Calculate dq/dp using the Implicit Function Theorem.
    # dq/dp = - (dF/dp) / (dF/dq)
    dF_dp_val = dF_dp(q_star, p_val)
    dF_dq_val = dF_dq(q_star, p_val)
    dq_dp = -dF_dp_val / dF_dq_val

    # Step 5: Calculate the elasticity of q with respect to p.
    # Elasticity ε = (dq/dp) * (p / q)
    elasticity = (dq_dp) * (p_val / q_star)

    # Step 6: Print the components of the elasticity calculation and the final result.
    print("The final elasticity is calculated as ε = (dq/dp) * (p / q*)")
    print(f"Optimal search intensity q* = {q_star:.5f}")
    print(f"Derivative dq/dp = {dq_dp:.5f}")
    print(f"Probability p = {p_val:.3f}")
    print(f"ε = {dq_dp:.3f} * ({p_val:.3f} / {q_star:.3f})")
    print(f"Calculated Elasticity: {elasticity:.3f}")

solve_and_calculate_elasticity()
<<< -0.341 >>>