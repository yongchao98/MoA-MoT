import numpy as np
from scipy.optimize import brentq

def solve_and_calculate_elasticity():
    """
    This function solves the microeconomic model to find the elasticity of optimal
    job search intensity (q) with respect to the probability of unemployment (p).
    """
    # Initial probability of unemployment
    p_initial = 0.2

    # The government's problem of maximizing the worker's expected utility
    # leads to a first-order condition F(q, p) = 0, which defines the
    # optimal search intensity q as an implicit function of p.
    # F(q, p) = exp(-2q) * [3 - 2q - p - 2pq + 2pq^2] - 2 + 2q = 0
    def F(q, p):
        """The implicit function F(q, p) = 0 from the FOC."""
        term1 = np.exp(-2 * q)
        term2 = 3 - 2*q - p - 2*p*q + 2*p*q**2
        return term1 * term2 - 2 + 2*q

    # Solve for the optimal q at p = 0.2 by finding the root of F(q, 0.2) = 0.
    # The search intensity q must be between 0 and 1.
    try:
        q_optimal = brentq(lambda q: F(q, p_initial), 0.0, 1.0)
    except ValueError:
        print("Error: Could not find the optimal q in the interval (0, 1).")
        return

    # To find dq/dp, we use the formula from the Implicit Function Theorem:
    # dq/dp = - (∂F/∂p) / (∂F/∂q)

    # ∂F/∂p = exp(-2q) * (-1 - 2q + 2q^2)
    def dF_dp(q, p):
        """Partial derivative of F with respect to p."""
        return np.exp(-2 * q) * (-1 - 2*q + 2*q**2)

    # ∂F/∂q = exp(-2q) * [-8 + 4q + 8pq - 4pq^2] + 2
    def dF_dq(q, p):
        """Partial derivative of F with respect to q."""
        term1 = np.exp(-2 * q)
        term2 = -8 + 4*q + 8*p*q - 4*p*q**2
        return term1 * term2 + 2

    # Evaluate the partial derivatives at the optimal q and initial p
    dF_dp_val = dF_dp(q_optimal, p_initial)
    dF_dq_val = dF_dq(q_optimal, p_initial)

    # Calculate the derivative of q with respect to p
    dq_dp = -dF_dp_val / dF_dq_val

    # Calculate the Elasticity: ε = (dq/dp) * (p / q)
    elasticity = dq_dp * (p_initial / q_optimal)

    # Output the results as requested
    print("The elasticity of optimal q with respect to p is calculated using the equation:")
    print("ε = (dq/dp) * (p / q)")
    print("\nFirst, we solve for the numbers in the equation:")
    print(f"The given probability of unemployment is p = {p_initial:.3f}")
    print(f"The calculated optimal search intensity is q = {q_optimal:.3f}")
    print(f"The calculated derivative of q with respect to p is dq/dp = {dq_dp:.3f}")
    print("\nPlugging these numbers into the elasticity equation:")
    print(f"ε = {dq_dp:.3f} * ({p_initial:.3f} / {q_optimal:.3f})")
    print(f"The final elasticity is {elasticity:.3f}")
    print(f"<<<{elasticity:.3f}>>>")

# Execute the function to print the solution
solve_and_calculate_elasticity()