import numpy as np

def solve_and_find_cost():
    """
    Solves the optimization problem to find the minimal cost per bit.
    """
    # The optimization problem reduces to finding the real root of p^3 + p - 1 = 0.
    # We define the polynomial coefficients for numpy.
    coeffs = [1, 0, 1, -1]

    # Find the roots of the polynomial.
    roots = np.roots(coeffs)

    # The optimal probability p* is the single real root of the equation.
    p_star = roots[np.isreal(roots)].real[0]

    # The minimal cost per bit of information is given by the formula C = -1 / log2(p*).
    # We use this formula for the final calculation.
    min_cost_per_bit = -1 / np.log2(p_star)

    # The problem asks for the constant factor in the total cost formula, which is this minimal cost per bit.
    # The final equation for the constant is C = -1 / log2(p*). We print the numbers involved.
    print(f"The equation to solve for the optimal probability is p^3 + p - 1 = 0.")
    print(f"The optimal probability is the real root, p* = {p_star:.6f}")
    print("\nThe final equation for the minimal cost per bit is: C = -1 / log2(p*)")
    print(f"Substituting the value of p*: C = -1 / log2({p_star:.6f})")
    print(f"\nThe minimal cost per bit of information is: {min_cost_per_bit:.3f}")

solve_and_find_cost()