import math
import numpy as np

def solve_and_print():
    """
    This function solves the problem by finding the constant 'c' which represents
    the minimal cost per bit of information.
    """
    
    # The problem reduces to solving the cubic equation x^3 + x - 1 = 0
    # where x = 2^(-1/c), and 'c' is the cost constant.
    # The coefficients of the polynomial x^3 + 0*x^2 + 1*x - 1 are [1, 0, 1, -1].
    coeffs = [1, 0, 1, -1]
    
    # We find the roots of this polynomial.
    roots = np.roots(coeffs)
    
    # We are interested in the single real root.
    real_root = None
    for root in roots:
        if np.isreal(root):
            real_root = np.real(root)
            break

    # Let's print the equation with the found root to show it holds.
    # We check if x^3 + x - 1 is indeed close to 0.
    equation_check_value = real_root**3 + real_root - 1
    
    print("The characteristic equation to solve is x^3 + x - 1 = 0.")
    print(f"The coefficients of the polynomial are: {coeffs[0]}, {coeffs[2]}, {coeffs[3]}")
    print(f"The real root of the equation is x = {real_root:.6f}")
    print(f"Checking the equation with the found root: {real_root:.6f}**3 + {real_root:.6f} - 1 = {equation_check_value:.6f}")
    print("-" * 20)
    
    # Now we calculate 'c' from x using the relation c = -ln(2)/ln(x).
    ln2 = math.log(2)
    ln_x = math.log(real_root)
    c = -ln2 / ln_x
    
    # Printing the final equation for 'c' with its components.
    print("The cost constant 'c' is calculated using the formula: c = -ln(2) / ln(x)")
    print(f"The value of ln(2) is: {ln2:.6f}")
    print(f"The value of ln(x) is: {ln_x:.6f}")
    print(f"The final calculation is: c = -({ln2:.6f}) / ({ln_x:.6f})")
    print("-" * 20)

    # The minimal number of coins per bit of information.
    print("The minimal cost per bit of information required is:")
    print(f"{c:.3f}")


solve_and_print()
