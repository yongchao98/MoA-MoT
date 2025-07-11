import numpy as np
from scipy.optimize import brentq

def solve_elasticity():
    """
    This function solves the economic problem to find the elasticity of optimal 
    job search intensity (q) with respect to the probability of unemployment (p).
    """
    
    # Given parameter
    p = 0.2

    # 1. Define the implicit function F(q, p) = 0 that characterizes the optimal q.
    # This equation is derived from the first-order conditions of the government's 
    # utility maximization problem subject to its budget constraint and the 
    # worker's incentive compatibility constraint.
    def F(q, p_val):
        """Implicit function F(q, p) = 0"""
        term1 = np.exp(-2 * q)
        term2 = p_val * (1 - q)**2 + (1 - p_val) * (1.5 - q)
        term3 = 1 - q
        return term1 * term2 - term3

    # 2. Solve for the optimal q (q_star) when p = 0.2
    # We need to find the root of F(q, 0.2) = 0. We know 0 < q < 1.
    # We use the brentq numerical solver for precision.
    try:
        q_star = brentq(lambda q: F(q, p), 1e-9, 1 - 1e-9)
    except ValueError:
        print("Could not find a root for q in the interval (0, 1).")
        return

    # 3. Use the Implicit Function Theorem to find dq/dp.
    # dq/dp = - (dF/dp) / (dF/dq)
    
    # Calculate the partial derivative dF/dp
    dF_dp = np.exp(-2 * q_star) * ((1 - q_star)**2 - (1.5 - q_star))

    # Calculate the partial derivative dF/dq
    # dF/dq = d/dq { exp(-2q) * [p(1-q)^2 + (1-p)(1.5-q)] - (1-q) }
    # Using the product rule and simplifying gives:
    exp_term = np.exp(-2 * q_star)
    h_q_p = p * (1 - q_star)**2 + (1 - p) * (1.5 - q_star)
    h_prime_q = -2 * p * (1 - q_star) - (1 - p)
    
    dF_dq = -2 * exp_term * h_q_p + exp_term * h_prime_q + 1

    # Calculate dq/dp
    dq_dp = -dF_dp / dF_dq
    
    # 4. Calculate the elasticity E = (dq/dp) * (p/q)
    elasticity = dq_dp * (p / q_star)
    
    print("This program calculates the elasticity of optimal job search intensity (q) with respect to the probability of unemployment (p).")
    print("\nThe final elasticity is calculated using the formula: Elasticity = (dq/dp) * (p / q*)")
    print("\nCalculated values for each component:")
    print(f"Optimal search intensity q* = {q_star:.5f}")
    print(f"Probability of unemployment p = {p}")
    print(f"Derivative dq/dp = {dq_dp:.5f}")
    
    print(f"\nPlugging these into the formula:")
    print(f"Elasticity = {dq_dp:.5f} * ({p} / {q_star:.5f})")
    print(f"Final Answer: Elasticity = {elasticity:.3f}")
    
    return elasticity

# Run the calculation and store the final answer
final_answer = solve_elasticity()

# Output the final answer in the required format
# print(f"\n<<<{final_answer:.3f}>>>")
print(f"\n<<<{final_answer:.3f}>>>")