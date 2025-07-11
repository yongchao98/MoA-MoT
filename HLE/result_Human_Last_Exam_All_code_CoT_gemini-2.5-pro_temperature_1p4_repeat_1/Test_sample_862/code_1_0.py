import numpy as np
from scipy.optimize import brentq

def solve_for_C():
    """
    This function solves the problem by finding the root of the transcendental
    equation and then calculating the constant C.
    """
    # The transcendental equation to solve is tan(x) - x - pi/2 = 0.
    # We are looking for the smallest positive root, which lies in (0, pi/2).
    def equation(x):
        return np.tan(x) - x - np.pi / 2

    # We use the brentq method for robust root finding in a given interval.
    # The interval must be chosen carefully to avoid the singularity of tan(x) at pi/2.
    # We search in a slighly reduced interval, e.g., (0.01, pi/2 - 0.01).
    try:
        search_interval_start = 0.01
        search_interval_end = np.pi / 2 - 0.01
        
        # Find the root x0
        x0 = brentq(equation, search_interval_start, search_interval_end)

        # The constant C is given by cot^2(x0/2)
        x0_div_2 = x0 / 2
        cot_val = 1 / np.tan(x0_div_2)
        C = cot_val**2
        
        # Print the results as requested
        print("The problem is reduced to solving the equation: tan(x) = x + pi/2")
        print(f"The smallest positive root is x0 = {x0:.10f}")
        print("\nThe constant C is found using the formula: C = cot^2(x0 / 2)")
        print("The numbers in the final equation are:")
        print(f"x0 / 2 = {x0_div_2:.10f}")
        print(f"cot(x0 / 2) = {cot_val:.10f}")
        print(f"C = ({cot_val:.10f})^2 = {C:.10f}")
        
    except (ValueError, RuntimeError) as e:
        print(f"An error occurred during the numerical computation: {e}")

solve_for_C()