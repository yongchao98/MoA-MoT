import math
import numpy as np

def solve_and_calculate_cost():
    """
    This function solves the problem by:
    1. Finding the root of the equation p^3 + p - 1 = 0.
    2. Calculating the minimal cost per bit K = -1 / log2(p).
    """

    # The problem of finding the optimal strategy reduces to solving the cubic equation p^3 + p - 1 = 0 for its real root p in (0, 1).
    # We can define the polynomial coefficients for numpy's root finder.
    # The equation is 1*p^3 + 0*p^2 + 1*p - 1 = 0
    coeffs = [1, 0, 1, -1]
    
    # Find all roots of the polynomial
    roots = np.roots(coeffs)
    
    # We are interested in the real root between 0 and 1.
    p_root = 0
    for r in roots:
        if np.isreal(r) and 0 < r < 1:
            p_root = np.real(r)
            break

    # The minimal cost per bit of information, K, is given by -1 / log2(p_root).
    # This value represents the leading constant factor of the total cost function Cost(n) / (n * log2(n)).
    cost_per_bit = -1 / math.log2(p_root)

    print("The optimal strategy is based on the solution to the equation p^3 + p - 1 = 0.")
    print(f"The relevant real root is p = {p_root:.6f}")
    
    # The final answer is the constant K, which is the cost per bit.
    # The "final equation" is the calculation of K from p.
    print("\nThe final cost calculation is:")
    log2_p = math.log2(p_root)
    print(f"K = -1 / log2(p)")
    print(f"K = {float(-1)} / log2({p_root:.6f})")
    print(f"K = {float(-1)} / {log2_p:.6f}")
    print(f"K = {cost_per_bit:.6f}")
    
    # Print the final answer rounded to 3 decimal places as requested.
    print(f"\nThe minimal cost per bit of information, rounded to 3 decimal places, is {cost_per_bit:.3f}.")

solve_and_calculate_cost()