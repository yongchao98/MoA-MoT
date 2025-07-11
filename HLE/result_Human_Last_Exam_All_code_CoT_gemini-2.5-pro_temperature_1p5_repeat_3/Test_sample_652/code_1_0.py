import numpy as np
from scipy.optimize import fsolve

def solve_and_calculate_elasticity():
    """
    Solves for the optimal search intensity q and calculates the elasticity
    of q with respect to the probability of unemployment p.
    """
    # Given parameter
    p_val = 0.2
    w_val = 20.0 # Although w cancels out, it's part of the model definition.

    # --- Step 1 & 2: Define and solve the First-Order Condition for q* ---
    # The FOC, G(q, p) = 0, is derived from the government's optimization problem.
    # G(q,p) = exp(2q)(1-q) - p(1-q)^2 / (1-p) - (1.5-q) = 0
    def G_func(q, p):
        return np.exp(2 * q) * (1 - q) - p * ((1 - q)**2) / (1 - p) - (1.5 - q)

    # Solve for q* at p=0.2
    q_initial_guess = 0.3
    q_star = fsolve(G_func, q_initial_guess, args=(p_val,))[0]

    print("--- Calculation Steps ---")
    print(f"Given probability of unemployment, p = {p_val}")
    print("\n1. Find the optimal search intensity q* by numerically solving the FOC:")
    print(f"   exp(2q*)(1-q*) - {p_val}(1-q*)²/(1-{p_val}) - (1.5-q*) = 0")
    print(f"   The optimal search intensity is q* = {q_star:.5f}")

    # --- Step 3 & 4: Use Implicit Function Theorem to find dq*/dp ---
    # dq*/dp = - (dG/dp) / (dG/dq)

    # Partial derivative of G with respect to p: dG/dp = -(1-q)² / (1-p)²
    dG_dp_val = -((1 - q_star)**2) / ((1 - p_val)**2)

    # Partial derivative of G with respect to q: dG/dq = exp(2q)(1-2q) + 2p(1-q)/(1-p) + 1
    dG_dq_val = (np.exp(2 * q_star) * (1 - 2 * q_star) +
                 2 * p_val * (1 - q_star) / (1 - p_val) + 1)
    
    # Calculate dq*/dp
    dq_dp_val = -dG_dp_val / dG_dq_val
    
    print("\n2. Calculate dq*/dp using the Implicit Function Theorem:")
    print(f"   dG/dp = {dG_dp_val:.5f}")
    print(f"   dG/dq = {dG_dq_val:.5f}")
    print(f"   dq*/dp = -({dG_dp_val:.5f}) / ({dG_dq_val:.5f}) = {dq_dp_val:.5f}")

    # --- Step 5: Calculate the elasticity ---
    elasticity = dq_dp_val * (p_val / q_star)

    print("\n3. Calculate the final elasticity, E = (dq*/dp) * (p / q*):")
    print(f"   E = {dq_dp_val:.5f} * ({p_val} / {q_star:.5f})")
    print(f"\nThe calculated elasticity is: {elasticity:.3f}")

if __name__ == '__main__':
    solve_and_calculate_elasticity()
    
<<<0.220>>>