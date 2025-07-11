import numpy as np

def solve_and_print_cost_constant():
    """
    This function solves for the optimal probability p, calculates the
    asymptotic cost constant C, and prints the results.
    """
    # The problem reduces to finding the constant C in the asymptotic cost formula:
    # Cost(n) ≈ C * (n * ln(n) - n)
    # C is derived from the root of the equation p^3 + p - 1 = 0.

    # Define the coefficients for the polynomial p^3 + 0*p^2 + 1*p - 1 = 0
    coeffs = [1, 0, 1, -1]

    # Find the roots of the polynomial
    roots = np.roots(coeffs)

    # The optimal probability p0 is the single real root of the equation.
    p0 = None
    for r in roots:
        if np.isreal(r):
            p0 = np.real(r)
            break

    if p0 is None:
        print("Could not find the real root.")
        return

    # The constant C is given by the formula -1 / ln(p0)
    constant_C = -1 / np.log(p0)

    print("The minimal cost is described by the asymptotic equation:")
    # The final equation for cost has two main numbers: C and -C
    # We output them as requested.
    print(f"Cost(n) ≈ {constant_C:.3f} * n*ln(n) - {constant_C:.3f} * n")
    print("\nThis is based on the following steps:")
    print("1. Solving the equation: p^3 + p - 1 = 0")
    print(f"   The real root is p = {p0:.7f}")
    print("2. Calculating the constant C = -1 / ln(p)")
    print(f"   C = -1 / {np.log(p0):.7f} ≈ {constant_C:.3f}")
    print("\nThe final answer is the value of this constant C.")


solve_and_print_cost_constant()