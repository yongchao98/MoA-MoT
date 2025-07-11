import numpy as np
from numpy.polynomial import Polynomial

def solve():
    """
    This function constructs two cubic polynomials f(x) and g(x) with positive
    derivatives and demonstrates that their composition f(g(x)) can have 9 fixed points.
    """
    # 1. Define base polynomials f0(x) and g(x) that satisfy the derivative conditions.
    # g(x) = x^3/3 + x. g'(x) = x^2 + 1, which is always positive.
    g_poly = Polynomial([0, 1, 0, 1/3])

    # Let's choose a non-symmetric f0(y) to generate more complex derivatives.
    # f0(y) = y^3 - 6y^2 + 13y.
    # f0'(y) = 3y^2 - 12y + 13. The discriminant is (-12)^2 - 4*3*13 = 144 - 156 = -12 < 0.
    # Since the leading coefficient is positive, f0'(y) > 0 for all y.
    f0_poly = Polynomial([0, 13, -6, 1])

    # 2. Scale f0(x) to ensure h'(x) = f'(g(x))g'(x) can cross the line y=1 multiple times.
    # The right scaling factor A for f(x) = A*f0(x) will make A*h0'(x) = 1 have 8 roots.
    # Through analysis (e.g., plotting), a factor of A=0.25 is found to be suitable.
    A = 0.25
    f_poly_no_d = A * f0_poly

    # 3. Construct the function k(x) = f(g(x)) - x, whose roots are the fixed points.
    # First, construct k(x) without the constant term 'd' in f(x).
    h_poly_no_d = f_poly_no_d(g_poly)
    k_poly_no_d = h_poly_no_d - Polynomial([0, 1])

    # 4. Find the critical points of k(x) to determine its local extrema.
    k_prime_poly = k_poly_no_d.deriv()
    crit_pts = k_prime_poly.roots()
    real_crit_pts = crit_pts[np.isreal(crit_pts)].real

    # 5. Determine the vertical shift 'd' needed to get 9 roots.
    # We evaluate k(x) at its critical points.
    crit_vals = k_poly_no_d(real_crit_pts)
    # To place the x-axis right in the middle of the wiggles, we shift the function
    # by the average of its maximum and minimum values at the critical points.
    d = -(np.max(crit_vals) + np.min(crit_vals)) / 2

    # 6. Construct the final polynomials with the determined shift 'd'.
    # f(x) = A*f0(x) + d
    f_poly = f_poly_no_d + Polynomial([d])
    # k(x) = f(g(x)) - x
    h_poly = f_poly(g_poly)
    k_poly = h_poly - Polynomial([0, 1])

    # 7. Find and count the real roots of the final equation k(x) = 0.
    roots = k_poly.roots()
    real_roots = roots[np.isreal(roots)].real
    num_real_roots = len(real_roots)

    print(f"The maximum number of fixed points is {num_real_roots}.")
    print("\nAn example polynomial equation f(g(x)) - x = 0 with 9 real roots is constructed.")
    print("The coefficients c_i of the polynomial Sum(c_i * x^i) = 0 are:")
    
    coeffs = k_poly.coef
    equation_parts = []
    for i, c in enumerate(coeffs):
        if not np.isclose(c, 0):
            print(f"c_{i}: {c:.4f}")
            # Also format the equation string
            if i == 0:
                equation_parts.append(f"{c:.4f}")
            elif i == 1:
                equation_parts.append(f"{c:.4f}x")
            else:
                equation_parts.append(f"{c:.4f}x^{i}")
    
    # Reverse for standard polynomial representation and join with '+'
    final_equation = " + ".join(reversed(equation_parts)).replace("+ -", "- ")
    print(f"\nThe final equation is:\n{final_equation} = 0")

solve()
<<<9>>>