import numpy as np
from scipy.optimize import brentq

def solve_and_calculate_elasticity():
    """
    Solves the model for the optimal q and calculates its elasticity with respect to p.
    """
    # Model parameters
    w = 20.0
    p_val = 0.2
    h = 1e-6 # A small change in p for numerical differentiation

    def get_q_from_b(b):
        """Calculates search intensity q given benefit b."""
        # Derived from the worker's FOC for search effort.
        # This relationship must hold: 0 < q < 1, which implies 20/e^2 < b < 20.
        if b <= 0 or b >= w:
            return -np.inf # Invalid q, return something to indicate failure
        return (np.log(w) - np.log(b)) / 2

    def government_foc(b, p, w_val):
        """
        Represents the government's First-Order Condition (FOC) for choosing b, d(E[U])/db.
        We find the root of this function to get the optimal b.
        """
        # Given b, calculate the worker's optimal q
        q = get_q_from_b(b)
        
        # Check for valid q
        if not (0 < q < 1):
             return np.inf

        # Government budget constraint: t = p*(1-q)*b
        tax = p * b * (1 - q)

        # Consumption when employed: c_e = w - t
        # Ensure consumption is positive for log utility
        if w_val - tax <= 0:
            return np.inf # Return a large number to guide the solver away

        # The FOC, derived by setting d(EU)/db = 0, is:
        # (1-q)/b - (1-p)*(1.5-q) / (w - t) = 0
        lhs = (1 - q) / b
        rhs = (1 - p) * (1.5 - q) / (w_val - tax)

        return lhs - rhs

    def get_optimal_q(p, w_val):
        """
        Solves the full model for a given p to find the optimal q.
        """
        # Search for optimal b in its valid range. Analysis shows the relevant root
        # is in the upper part of the range for p=0.2.
        b_lower_bound = 10.0
        b_upper_bound = w_val - 1e-5  # Stay just inside the b < w domain

        try:
            # Find the optimal b* that makes the FOC zero
            optimal_b = brentq(government_foc, b_lower_bound, b_upper_bound, args=(p, w_val))
        except ValueError:
            # This can happen if the FOC does not cross zero in the interval,
            # indicating a potential issue or a corner solution under different parameters.
            print("Error: Could not find an optimal b in the specified interval.")
            return None

        # Calculate the resulting optimal q* from the optimal b*
        optimal_q = get_q_from_b(optimal_b)
        return optimal_q

    # 1. Calculate optimal q at p
    q_at_p = get_optimal_q(p_val, w)
    if q_at_p is None:
        return

    # 2. Calculate optimal q at p+h
    q_at_p_plus_h = get_optimal_q(p_val + h, w)
    if q_at_p_plus_h is None:
        return

    # 3. Calculate the numerical derivative dq/dp
    derivative_q_p = (q_at_p_plus_h - q_at_p) / h

    # 4. Calculate the elasticity
    elasticity = derivative_q_p * (p_val / q_at_p)

    print("--- Calculation of the Elasticity ---")
    print(f"Optimal search intensity q* at p = {p_val:.3f} is: {q_at_p:.3f}")
    print(f"The derivative dq*/dp at p = {p_val:.3f} is approximately: {derivative_q_p:.3f}")
    print("\nElasticity of q* with respect to p is given by the formula: (dq*/dp) * (p / q*)")
    print("Plugging in the numbers:")
    print(f"Elasticity = {derivative_q_p:.3f} * ({p_val:.3f} / {q_at_p:.3f})")
    print(f"\nThe calculated elasticity is: {elasticity:.3f}")

solve_and_calculate_elasticity()