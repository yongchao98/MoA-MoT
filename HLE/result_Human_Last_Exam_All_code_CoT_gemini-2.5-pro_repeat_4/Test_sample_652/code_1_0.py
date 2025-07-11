import numpy as np
from scipy.optimize import fsolve

def solve_and_calculate_elasticity():
    """
    Solves the microeconomic model to find the elasticity of optimal job search
    intensity (q) with respect to the probability of unemployment (p).
    """
    # Given parameters
    p = 0.2
    w = 20

    # The government's optimization problem can be reduced to a single implicit
    # equation that defines the optimal search intensity q* as a function of p.
    # This equation, which we'll call H(q, p), is derived from the first-order
    # conditions of the government's problem.
    # H(q, p) = (1-q)exp(2q) - (p/(1-p))(1-q)^2 - 1.5 + q = 0
    def H(q, p_val):
        return (1 - q) * np.exp(2 * q) - (p_val / (1 - p_val)) * (1 - q)**2 - 1.5 + q

    # --- Step 1: Solve for the optimal q (q_star) at p = 0.2 ---
    # We need a function of q only for the numerical solver.
    func_q_for_solver = lambda q: H(q, p)
    # An initial guess for q is needed for the solver. q is a probability between 0 and 1.
    q_initial_guess = 0.3
    # Use fsolve to find the root of the equation H(q, 0.2) = 0.
    q_star = fsolve(func_q_for_solver, q_initial_guess)[0]

    # --- Step 2: Calculate the derivative dq/dp using the Implicit Function Theorem ---
    # The theorem states dq/dp = - (∂H/∂p) / (∂H/∂q).
    # We define functions for these partial derivatives.
    def dH_dp(q, p_val):
        """Partial derivative of H with respect to p."""
        return - (1 / (1 - p_val)**2) * (1 - q)**2

    def dH_dq(q, p_val):
        """Partial derivative of H with respect to q."""
        term1 = (1 - 2*q) * np.exp(2 * q)  # Derivative of (1-q)exp(2q)
        term2 = (2 * p_val * (1 - q)) / (1 - p_val)  # Derivative of -(p/(1-p))(1-q)^2
        term3 = 1  # Derivative of q
        return term1 + term2 + term3

    # Calculate the numerical values of the partial derivatives at (q_star, p).
    dH_dp_val = dH_dp(q_star, p)
    dH_dq_val = dH_dq(q_star, p)

    # Calculate dq/dp.
    dq_dp = -dH_dp_val / dH_dq_val

    # --- Step 3: Calculate the elasticity of q with respect to p ---
    elasticity = dq_dp * (p / q_star)

    # --- Output the results step-by-step ---
    print("The optimal search intensity q is implicitly defined by the equation:")
    print("H(q, p) = (1-q)exp(2q) - [p/(1-p)](1-q)^2 - 1.5 + q = 0\n")

    print(f"Step 1: Solve for the optimal q (q*) when p = {p}")
    print(f"Solving numerically, we find q* = {q_star:.5f}\n")

    print("Step 2: Calculate the derivative dq/dp using the Implicit Function Theorem.")
    print("dq/dp = - (∂H/∂p) / (∂H/∂q)")
    print(f"The partial derivative ∂H/∂p at (q*, p) is: {dH_dp_val:.5f}")
    print(f"The partial derivative ∂H/∂q at (q*, p) is: {dH_dq_val:.5f}")
    print(f"Therefore, dq/dp = -({dH_dp_val:.5f}) / ({dH_dq_val:.5f}) = {dq_dp:.5f}\n")

    print("Step 3: Calculate the final elasticity of q with respect to p.")
    print("Elasticity = (dq/dp) * (p / q*)")
    print(f"Elasticity = {dq_dp:.5f} * ({p} / {q_star:.5f})")
    print(f"The calculated elasticity is: {elasticity:.3f}")

    # Final answer in the required format
    print(f"\n<<<{elasticity:.3f}>>>")

# Run the calculation and print the results
solve_and_calculate_elasticity()