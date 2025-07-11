import numpy as np
from scipy.optimize import brentq

def solve_and_calculate_elasticity():
    """
    This function solves the microeconomic model to find the elasticity of optimal
    job search intensity (q) with respect to the probability of unemployment (p).
    """
    # Given parameters
    p = 0.2
    w = 20

    # Step 1: Find the optimal q for p=0.2
    # The optimal q is the root of the equation H(q, p) = 0. For p=0.2, this simplifies to:
    # exp(-2q) * (q^2 - 6q + 7) - 4 + 4q = 0
    # We define this function to find its root.
    def f_q_for_p02(q):
        return np.exp(-2 * q) * (q**2 - 6 * q + 7) - 4 + 4 * q

    # Numerically solve for q in the interval (0, 1). Based on analysis, the root is in (0, 0.5).
    try:
        q_star = brentq(f_q_for_p02, 0, 0.5)
    except ValueError:
        print("Error: Could not find the root for q in the specified interval.")
        return

    # Step 2: Define the partial derivatives of the implicit function H(q, p).
    # H(q, p) = (1-p)*exp(2q)*(1-q) - (1-p)*(1.5-q) - p*(1-q)^2 = 0
    
    # Partial derivative of H with respect to p
    def dH_dp(q, p):
        return -np.exp(2 * q) * (1 - q) + 0.5 + q - q**2

    # Partial derivative of H with respect to q
    def dH_dq(q, p):
        return (1 - p) * np.exp(2 * q) * (1 - 2 * q) + (1 - p) + 2 * p * (1 - q)

    # Step 3: Calculate the numerical values of the partial derivatives at (q_star, p)
    dH_dp_val = dH_dp(q_star, p)
    dH_dq_val = dH_dq(q_star, p)
    
    # Step 4: Calculate dq/dp using the Implicit Function Theorem
    # dq/dp = - (dH/dp) / (dH/dq)
    if dH_dq_val == 0:
        print("Error: Division by zero. dH/dq is zero.")
        return
    dq_dp = -dH_dp_val / dH_dq_val

    # Step 5: Calculate the elasticity
    # ε = (dq/dp) * (p/q)
    if q_star == 0:
        print("Error: Division by zero. q_star is zero.")
        return
    elasticity = dq_dp * (p / q_star)

    # Output the results
    print(f"For p = {p}:")
    print(f"The optimal job search intensity is q* = {q_star:.5f}")
    print("\nCalculating the elasticity ε = (dq/dp) * (p/q*):")
    print(f"The derivative dq/dp is calculated to be {dq_dp:.5f}")
    
    # As requested, output the numbers in the final equation
    print(f"\nThe final elasticity calculation is:")
    print(f"ε = {dq_dp:.5f} * ({p:.1f} / {q_star:.5f})")
    
    print(f"\nThe elasticity of optimal q with respect to p is: {elasticity:.3f}")
    
    # Final answer in the specified format
    print(f"\n<<<{elasticity:.3f}>>>")

if __name__ == "__main__":
    solve_and_calculate_elasticity()
