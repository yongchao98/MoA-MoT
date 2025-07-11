import numpy as np
from scipy.optimize import brentq

def solve_elasticity():
    """
    This function solves the problem step-by-step as outlined above.
    1. Finds the optimal q for p=0.2.
    2. Calculates the partial derivatives needed for the implicit function theorem.
    3. Computes the elasticity of q with respect to p.
    """
    # The probability of becoming unemployed
    p = 0.2

    # The government's first-order condition H(q, p) = 0
    # This equation implicitly defines the optimal q* as a function of p.
    def H(q, p):
        return np.exp(-2 * q) * (3 - p - 2 * q * (1 + p) + 2 * p * q**2) - 2 * (1 - q)

    # We solve for q* at p=0.2. Analysis of H(q, 0.2) shows a root between
    # 0 and 0.5 that corresponds to a utility maximum.
    try:
        q_star = brentq(H, 0.0, 0.5, args=(p,))
    except ValueError:
        print("Could not find a root in the specified interval.")
        return

    # Partial derivative of H with respect to p
    def dH_dp(q, p):
        return np.exp(-2 * q) * (2 * q**2 - 2 * q - 1)

    # Partial derivative of H with respect to q
    def dH_dq(q, p):
        # The SOC is equivalent to dH/dq < 0 at the root.
        return 4 * q - 2 - 2 * np.exp(-2 * q) * (1 + p - 2 * p * q)

    # Evaluate the partial derivatives at (q_star, p)
    dH_dp_val = dH_dp(q_star, p)
    dH_dq_val = dH_dq(q_star, p)

    # Calculate dq/dp using the Implicit Function Theorem
    dq_dp = - (dH_dp_val / dH_dq_val)

    # Calculate the elasticity E = (dq/dp) * (p/q)
    elasticity = dq_dp * (p / q_star)

    # --- Output the results ---
    print(f"Step 1: Find the optimal search intensity q* for p = {p}")
    print(f"The optimal q* that maximizes utility is: {q_star:.4f}\n")
    
    print("Step 2: Calculate the components of the elasticity formula E = (dq/dp) * (p/q*)")
    print(f"The partial derivative dH/dp at (q*, p) is: {dH_dp_val:.4f}")
    print(f"The partial derivative dH/dq at (q*, p) is: {dH_dq_val:.4f}\n")
    
    print(f"Step 3: Calculate dq/dp using the Implicit Function Theorem: dq/dp = - (dH/dp / dH/dq)")
    print(f"dq/dp = -({dH_dp_val:.4f} / {dH_dq_val:.4f}) = {dq_dp:.4f}\n")
    
    print("Step 4: Calculate the final elasticity.")
    print("The final equation is: Elasticity = (dq/dp) * (p / q*)")
    print(f"Plugging in the numbers: Elasticity = {dq_dp:.4f} * ({p} / {q_star:.4f})")
    print(f"The elasticity of optimal q with respect to p is: {elasticity:.3f}")
    
    print(f"\n<<<{elasticity:.3f}>>>")

solve_elasticity()