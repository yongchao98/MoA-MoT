import numpy as np
import math

def solve_and_print_cost():
    """
    Solves for the optimal strategy and calculates the minimal cost coefficient for sorting.
    """
    # The problem of finding the minimal cost per bit of information boils down to
    # minimizing the function C(p) = (3 - 2p) / (-p*log2(p) - (1-p)*log2(1-p)).
    # The minimum is found when p is the real root of the polynomial p^3 + p - 1 = 0.

    # Coefficients for the polynomial p^3 + 0*p^2 + 1*p - 1 = 0
    coeffs = [1, 0, 1, -1]

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # The equation has one real root and two complex conjugate roots. We need the real root.
    p0 = roots[np.isreal(roots)].real[0]

    # The minimal cost per bit of information is C_min = -1 / log2(p0)
    # This is the coefficient for the leading term of the total cost function.
    c_min = -1 / math.log2(p0)

    # The total cost is asymptotically C_min * n * log2(n)
    print("The minimal number of coins to sort the array is given by the asymptotic formula:")
    # We output each number in the final equation as requested.
    print(f"Total Cost = {c_min:.3f} * n * log_2(n)")

solve_and_print_cost()

# For the final answer, we calculate the coefficient C_min with higher precision.
coeffs = [1, 0, 1, -1]
roots = np.roots(coeffs)
p0 = roots[np.isreal(roots)].real[0]
c_min = -1 / math.log2(p0)
# We are asked for the answer (the coefficient) up to 3 decimal places.
# print(f"<<<{c_min:.3f}>>>")