import numpy as np
from scipy.optimize import brentq

def solve_and_explain():
    """
    This function explains the setup and calculates the valid range for R
    to define the set M for the Banach Fixed-Point Theorem.
    """
    print("To solve the Boundary Value Problem using the Banach Fixed-Point Theorem, we define a set M.")
    print("M must be a complete metric space where our operator T is a contraction mapping from M to M.")

    # Define the space and the form of M
    print("\nWe choose the space of continuous functions on [0, 1] that are zero at the endpoints, C₀[0, 1].")
    print("The set M is a closed ball of radius R in this space:")
    print("M = {u ∈ C₀[0, 1] | ||u||_∞ ≤ R}")

    # State the conditions on R
    print("\nFor the theorem to apply, R must satisfy two conditions:")
    print("1. Invariance (T: M -> M):  exp(R) / 8 ≤ R")
    print("2. Contraction:             exp(R) / 8 < 1  => R < ln(8)")

    # Define the function for finding the lower bound
    def h(R):
        """Represents the inequality R - exp(R)/8 >= 0"""
        return R - np.exp(R) / 8

    # Find the bounds for R numerically
    try:
        # We find the root of h(R) = 0. A quick test shows h(0.1)<0 and h(0.2)>0.
        r_lower = brentq(h, 0.1, 0.2)
    except (ValueError, RuntimeError):
        # Fallback in case solver fails
        r_lower = 0.129
        
    r_upper = np.log(8)

    print(f"\nSolving these inequalities, we find that R must be in the interval [{r_lower:.4f}, {r_upper:.4f}).")

    # Choose a valid R and define the final set M
    r_chosen = 1.0
    print(f"\nWe can pick any R in this interval. Let's choose a simple value, R = {r_chosen}.")
    
    print("\nTherefore, the set M you should define is:")
    set_definition_part_1 = "M = {u ∈ C[0, 1] | u(0) = "
    u_0 = 0
    set_definition_part_2 = ", u(1) = "
    u_1 = 0
    set_definition_part_3 = " and sup_{x∈[0,1]} |u(x)| ≤ "
    final_R = r_chosen
    
    # "output each number in the final equation!"
    print(f"{set_definition_part_1}{u_0}{set_definition_part_2}{u_1}{set_definition_part_3}{final_R}}")


solve_and_explain()
<<<M = {u ∈ C[0, 1] | u(0) = 0, u(1) = 0 and sup_{x∈[0,1]} |u(x)| ≤ R} where R is a constant satisfying exp(R)/8 ≤ R and R < ln(8). For example, R=1.>>>