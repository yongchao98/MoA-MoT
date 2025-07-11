import numpy as np
from scipy.optimize import fsolve

def solve_alpha():
    """
    This function solves for the value of alpha based on the problem description.
    """
    # Step 1: Define constants from the problem statement
    T = np.log(10)
    B_val = 0.5 * 10**20 / 99**2

    # Step 2: Define intermediate constants K2 and K3 based on T
    # K2 = 1 - exp(-2T)
    # K3 = 1 - exp(-3T)
    K2 = 1 - np.exp(-2 * T)
    K3 = 1 - np.exp(-3 * T)

    # Step 3: Define the algebraic equation for alpha.
    # The derivation leads to the equation: C1 * alpha^5 - C2 * alpha^8 = B
    # where:
    # C1 = 12 / (K3 * K2^4)
    # C2 = 32 / K2^8
    # Let's define a function f(alpha) = C1 * alpha^5 - C2 * alpha^8 - B = 0

    C1 = 12 / (K3 * K2**4)
    C2 = 32 / (K2**8)

    # The final equation to solve for alpha is:
    # C1 * alpha^5 - C2 * alpha^8 - B = 0
    
    print("The final equation to be solved is of the form: C1 * alpha^5 - C2 * alpha^8 = B")
    print(f"where C1 = {C1}")
    print(f"      C2 = {C2}")
    print(f"      B = {B_val}")
    
    # We define the function to find its root
    def equation_for_alpha(alpha):
        # We only search for positive alpha as per the problem statement (alpha > 0)
        if alpha <= 0:
            return np.inf # Return a large number for non-physical values
        return C1 * alpha**5 - C2 * alpha**8 - B_val

    # Step 4: Solve the equation numerically.
    # A consistency check from the derivation shows that for real solutions to exist,
    # alpha^3 must be less than (3 * K2^4) / (8 * K3).
    # This implies alpha must be small (less than 1).
    # We will therefore look for a small positive root.
    initial_guess = 0.1
    alpha_solution = fsolve(equation_for_alpha, initial_guess)[0]

    # Output the final answer
    print("\nThe calculated value of alpha is:")
    print(alpha_solution)

solve_alpha()