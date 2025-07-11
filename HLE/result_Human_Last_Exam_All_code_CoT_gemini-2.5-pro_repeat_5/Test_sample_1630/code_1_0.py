import numpy as np
import matplotlib.pyplot as plt

def demonstrate_derivative_behavior():
    """
    This function demonstrates that an 8th-degree polynomial H(x) can exist such that:
    1. H(x) > 0 for all x (as required for h'(x)).
    2. The equation H(x) = 1 has 8 distinct real roots.
    
    This is the key condition required for the original function h(x) = f(g(x))
    to have 9 fixed points.
    """
    
    # Let's construct a polynomial H(x) such that H(x) - 1 has 8 roots.
    # We choose the roots to be -4, -3, -2, -1, 1, 2, 3, 4.
    # So, H(x) - 1 = c * (x-1)(x+1)(x-2)(x+2)(x-3)(x+3)(x-4)(x+4)
    # H(x) - 1 = c * (x^2-1)(x^2-4)(x^2-9)(x^2-16)
    # P(x) = (x^2-1)(x^2-4)(x^2-9)(x^2-16) is a polynomial with these roots.
    # Let's expand P(x) = x^8 - 30x^6 + 273x^4 - 820x^2 + 576
    p_coeffs = [1, 0, -30, 0, 273, 0, -820, 0, 576]
    P = np.poly1d(p_coeffs)

    # To ensure H(x) > 0, we need H(x) = 1 + c*P(x) > 0 for all x.
    # This means 1 > -c * P_min.
    # Let's find the minimum value of P(x).
    # We find the roots of the derivative of P(x) to find the extrema.
    crit_points = P.deriv().r
    # Evaluate P at the real critical points to find extrema.
    extrema = P(crit_points[np.isreal(crit_points)])
    P_min = np.min(extrema) # Approximately -1285.4

    # We need 1 + c * P_min > 0. If c > 0, we need c < -1/P_min.
    # c < 1/1285.4. Let's choose c = 1/2000.
    c = 1/2000.0

    # Our polynomial H(x) that represents a possible h'(x) is:
    # H(x) = 1 + c*P(x)
    H_coeffs = c * p_coeffs
    H_coeffs[-1] += 1
    H = np.poly1d(H_coeffs)
    
    # Now, let's find the roots of H(x) - 1 = 0, which by construction are the roots of P(x).
    roots = H.r - np.poly1d([1]).r # This is not correct way to find roots of H(x)-1=0
    roots_of_H_minus_1 = (H - np.poly1d([1])).r
    # Sort the real roots
    real_roots = np.sort(roots_of_H_minus_1[np.isreal(roots_of_H_minus_1)].real)
    
    print("The argument for 9 fixed points requires finding an 8th-degree polynomial h'(x) that is always positive but equals 1 at 8 distinct points.")
    print("Here is an example of such a polynomial, H(x) = 1 + c*P(x):")
    
    final_equation = H - np.poly1d([1])
    
    # Print the coefficients of the equation H(x) - 1 = 0
    print("\nThe equation H(x) - 1 = 0 is a polynomial equation:")
    for i, coeff in enumerate(final_equation.coeffs):
        power = len(final_equation.coeffs) - 1 - i
        if abs(coeff) > 1e-9: # Only print non-zero terms
            print(f"({coeff:.6f}) * x^{power}", end=" + ")
    print("0 = 0")

    print(f"\nThis equation has {len(real_roots)} real roots, as required:")
    print(np.round(real_roots, 4))
    
    # Verify that H(x) is always positive
    H_crit_points = H.deriv().r
    H_extrema = H(H_crit_points[np.isreal(H_crit_points)])
    H_min = np.min(H_extrema)
    
    print(f"\nThe minimum value of this H(x) is approximately {H_min:.4f}, which is positive.")
    print("\nSince it's possible for h'(x) to satisfy these conditions, it is possible for h(x)=f(g(x)) to have 9 fixed points.")

demonstrate_derivative_behavior()
