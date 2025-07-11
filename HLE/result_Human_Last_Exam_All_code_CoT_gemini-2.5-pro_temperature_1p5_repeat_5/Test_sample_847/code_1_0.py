import numpy as np
import math

def solve_for_minimal_cost_coefficient():
    """
    This function calculates the minimal cost coefficient for sorting the array.

    The problem reduces to finding the most cost-effective way to gain information.
    The cost per bit 'k' for Type 2 questions is minimized when p^3 + p - 1 = 0,
    where p is the probability of a 'yes' answer.
    The cost 'k' is then calculated as k = -1 / log2(p).
    """

    # Coefficients of the polynomial p^3 + 0*p^2 + 1*p - 1 = 0
    coeffs = [1, 0, 1, -1]

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # The real root is the probability 'p' we are interested in.
    p_y = roots[np.isreal(roots)].real[0]

    # Calculate the minimal cost per bit, k.
    k = -1 / math.log2(p_y)

    # As per the instruction to "output each number in the final equation",
    # we show the final calculation for k.
    print(f"The optimal strategy involves designing questions where the probability of a 'yes' answer is the real root of p^3 + p - 1 = 0.")
    print(f"The real root is p = {p_y:.6f}")
    print(f"The minimal cost per bit is k = -1 / log2(p).")
    print(f"Substituting p, the equation is: k = -1 / log2({p_y:.6f}) = {k:.6f}")
    
    print("\nSince this cost per bit (~1.813) is less than the 2 coins/bit from comparison questions, this is the optimal rate.")
    print(f"The minimal coefficient, rounded to 3 decimal places, is: {k:.3f}")


solve_for_minimal_cost_coefficient()