import numpy as np
from scipy.optimize import root_scalar

def solve_elasticity():
    """
    This function calculates the elasticity of the optimal re-employment probability q
    with respect to the unemployment probability p.
    """
    # Define the problem parameter
    p = 0.2

    # Define the implicit function F(q, p) = 0 that determines the optimal q.
    # F(q, p) = 1.5 - q + (p/(1-p)) * (1-q)**2 - (1-q) * exp(2q)
    def F(q, p_val):
        return 1.5 - q + (p_val / (1 - p_val)) * (1 - q)**2 - (1 - q) * np.exp(2 * q)

    # Use a numerical solver to find the optimal q* for p = 0.2.
    # We search for a root in the valid probability interval (0, 1).
    try:
        solution = root_scalar(F, args=(p,), bracket=[1e-6, 1 - 1e-6])
        q_star = solution.root
    except ValueError:
        print("Error: Could not find a root in the specified interval.")
        return

    # Define the partial derivative of F with respect to q.
    # dF/dq = -1 - 2*(p/(1-p))*(1-q) - exp(2q)*(1 - 2q)
    def dF_dq(q, p_val):
        return -1 - 2 * (p_val / (1 - p_val)) * (1 - q) - np.exp(2 * q) * (1 - 2 * q)

    # Define the partial derivative of F with respect to p.
    # dF/dp = (1 / (1-p)**2) * (1-q)**2
    def dF_dp(q, p_val):
        return (1 / (1 - p_val)**2) * (1 - q)**2

    # Calculate the values of the partial derivatives at (q*, p).
    dFdq_val = dF_dq(q_star, p)
    dFdp_val = dF_dp(q_star, p)

    # Calculate dq/dp using the Implicit Function Theorem.
    dq_dp = -dFdp_val / dFdq_val

    # Calculate the elasticity of q with respect to p.
    elasticity = dq_dp * (p / q_star)

    print(f"The elasticity of optimal q with respect to p is {elasticity:.3f}.")

solve_elasticity()