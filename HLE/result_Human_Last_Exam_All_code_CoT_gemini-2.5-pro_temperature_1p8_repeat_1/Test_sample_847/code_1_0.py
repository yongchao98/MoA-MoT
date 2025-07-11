import numpy as np
import math

def solve():
    """
    This function calculates the minimal cost per bit of information for the given problem.

    The cost is determined by finding the optimal strategy between two types of questions.
    This leads to solving a polynomial equation y^3 - y^2 - 2y - 1 = 0 to find the optimal
    split ratio y0 for general questions. The minimal cost per bit, C, is then derived
    from this root using the formula C = 2 / log2(y0).
    """

    # Coefficients of the polynomial y^3 - y^2 - 2y - 1 = 0
    # The coefficients are for y^3, y^2, y^1, y^0
    coeffs = [1, -1, -2, -1]

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # We are interested in the single positive real root
    y0 = None
    for r in roots:
        # Check if the root is real (imaginary part is close to zero)
        if np.isreal(r):
            # Check if the real part is positive
            if r.real > 0:
                y0 = r.real
                break

    if y0 is None:
        print("Could not find the positive real root.")
        return

    # Calculate C = 2 / log2(y0)
    log2_y0 = math.log2(y0)
    C = 2 / log2_y0

    # Output the steps of the final calculation
    print("The optimal questioning strategy leads to solving the characteristic equation:")
    print("y^3 - y^2 - 2y - 1 = 0")
    print("\nThe positive real root of this equation determines the optimal split ratio:")
    print(f"y0 = {y0:.3f}")
    
    print("\nThe minimal cost per bit of information, C, is calculated as:")
    print("C = 2 / log2(y0)")
    
    print("\nSubstituting the values:")
    print(f"C = 2 / log2({y0:.3f})")
    print(f"C = 2 / {log2_y0:.3f}")
    print(f"C = {C:.3f}")

solve()
