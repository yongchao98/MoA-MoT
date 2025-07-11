import numpy as np
from scipy.optimize import brentq

def solve_and_calculate_elasticity():
    """
    This function solves for the optimal search intensity q* and calculates its
    elasticity with respect to the unemployment probability p.
    """
    # Define the problem parameter p
    p = 0.2

    # Step 1: Define the First-Order Condition (FOC) for the optimal q*.
    # The government's problem simplifies to maximizing EU with respect to q,
    # and the FOC is F(q,p) = 0.
    # F(q,p) = [(1-p)(3-2q) + 2p(1-q)^2] * exp(-2q) - 2(1-q) = 0
    def foc_equation(q, p_val):
        term1 = (1 - p_val) * (3 - 2 * q)
        term2 = 2 * p_val * (1 - q)**2
        return (term1 + term2) * np.exp(-2 * q) - 2 * (1 - q)

    # Step 2: Solve for q* at p = 0.2 using a numerical root finder.
    # q is a probability, so it must be in the interval (0, 1).
    try:
        q_star = brentq(lambda q: foc_equation(q, p), 1e-5, 1 - 1e-5)
    except ValueError:
        print("Could not find a unique solution for q* in the interval (0,1).")
        return

    # Step 3: Calculate the derivative dq*/dp using the Implicit Function Theorem.
    # dq*/dp = - (dF/dp) / (dF/dq), where F is the FOC.
    
    # Define the partial derivative of F with respect to p:
    # dF/dp = (2q^2 - 2q - 1) * exp(-2q)
    def dF_dp(q):
        return (2 * q**2 - 2 * q - 1) * np.exp(-2 * q)

    # Define the partial derivative of F with respect to q:
    # dF/dq = exp(-2q) * [-4pq^2 + 8pq + 4q - 8] + 2
    def dF_dq(q, p_val):
        bracket_term = -4 * p_val * q**2 + 8 * p_val * q + 4 * q - 8
        return np.exp(-2 * q) * bracket_term + 2

    # Calculate the numerical values of the partial derivatives at (q*, p).
    dF_dp_val = dF_dp(q_star)
    dF_dq_val = dF_dq(q_star, p)

    # Calculate dq*/dp
    dq_dp_val = -dF_dp_val / dF_dq_val

    # Step 4: Calculate the elasticity E = (dq*/dp) * (p / q*).
    elasticity = dq_dp_val * (p / q_star)

    # Print the final result, including the numbers used in the elasticity formula.
    print("The elasticity is calculated using the formula: E = (dq*/dp) * (p / q*)")
    print(f"Optimal search intensity q* = {q_star:.5f}")
    print(f"The derivative dq*/dp = {dq_dp_val:.5f}")
    print(f"The unemployment probability p = {p:.1f}")
    print(f"\nThe elasticity calculation is: E = {dq_dp_val:.5f} * ({p:.1f} / {q_star:.5f})")
    print(f"The resulting elasticity, rounded to three decimals, is: {elasticity:.3f}")


solve_and_calculate_elasticity()
<<< -0.329 >>>