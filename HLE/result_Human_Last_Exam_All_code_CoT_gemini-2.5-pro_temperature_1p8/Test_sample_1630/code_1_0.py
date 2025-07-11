import numpy as np
from numpy.polynomial import polynomial as P

def solve():
    """
    This function demonstrates that a strictly increasing polynomial of degree 9
    can have 9 fixed points.
    """
    # 1. Define a polynomial P9(x) with 9 distinct real roots: -4, -3, -2, -1, 0, 1, 2, 3, 4.
    # The polynomial is P9(x) = x * (x^2-1) * (x^2-4) * (x^2-9) * (x^2-16)
    roots_P9 = np.array([-4, -3, -2, -1, 0, 1, 2, 3, 4])
    # The coefficients can be obtained from the roots.
    # The numpy.polynomial.Polynomial class handles this conversion.
    poly_P9 = P.Polynomial.fromroots(roots_P9)

    # 2. Compute the derivative P9'(x).
    poly_P9_prime = poly_P9.deriv()

    # 3. Find the minimum value of P9'(x).
    # The minima of P9'(x) occur at the roots of its derivative, P9''(x).
    poly_P9_double_prime = poly_P9_prime.deriv()
    crit_pts_of_prime = poly_P9_double_prime.roots()
    # We are interested in real roots only.
    real_crit_pts_of_prime = crit_pts_of_prime[np.isreal(crit_pts_of_prime)].real
    # Evaluate P9'(x) at these critical points to find its local extrema.
    values_at_extrema = poly_P9_prime(real_crit_pts_of_prime)
    min_val_P9_prime = np.min(values_at_extrema)

    # 4. Choose a constant C to make H'(x) > 0.
    # We need C * P9'(x) + 1 > 0.
    # Let's choose C = -1 / (2 * min_val_P9_prime). Since min_val is negative, C is positive.
    C = -1.0 / (2.0 * min_val_P9_prime)

    # 5. Define the polynomial H(x) = C * P9(x) + x.
    # H'(x) = C * P9'(x) + 1. Its minimum value will be C*min_val + 1 = 0.5 > 0.
    # Thus H(x) is strictly increasing.
    # In polynomial form:
    poly_x = P.Polynomial([0, 1])
    poly_H = C * poly_P9 + poly_x

    # 6. Find the fixed points of H(x).
    # This means solving H(x) = x, which is C * P9(x) = 0.
    # The roots are the same as the roots of P9(x).
    fixed_points = poly_H.roots() - poly_x.roots() # This is conceptually what happens. Numerically solve C*P9(x) = 0
    final_poly_to_solve = C * poly_P9
    fixed_points = final_poly_to_solve.roots()

    # Filter for real roots and sort them
    real_fixed_points = sorted(np.round(fixed_points[np.isreal(fixed_points)].real, 5))
    
    print("A strictly increasing polynomial H(x) of degree 9 was constructed.")
    min_H_prime = C * min_val_P9_prime + 1
    print(f"The minimum value of its derivative H'(x) is {min_H_prime:.2f}, which is positive.")
    print("The fixed points of H(x) are the solutions to H(x) = x.")
    print("These are the roots of the equation: C * (x^9 - 20x^7 + ...) = 0")
    print("The 9 fixed points found are:")
    for point in real_fixed_points:
        print(point)
    print(f"\nThe maximum number of fixed points is {len(real_fixed_points)}.")

solve()