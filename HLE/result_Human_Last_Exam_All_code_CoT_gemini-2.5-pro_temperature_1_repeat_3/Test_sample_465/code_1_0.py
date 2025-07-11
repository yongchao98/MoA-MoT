import numpy as np

def solve_and_check_roots():
    """
    This function demonstrates that chromatic roots can be non-real (Statement B).
    It calculates the roots of the chromatic polynomial for the 4-cycle graph (C4).
    The chromatic polynomial for C4 is P(x) = x^4 - 4x^3 + 6x^2 - 3x.
    """
    # Coefficients of the polynomial x^4 - 4x^3 + 6x^2 - 3x + 0
    # in descending order of power.
    coeffs = [1, -4, 6, -3, 0]

    # Calculate the roots of the polynomial
    roots = np.roots(coeffs)

    print("The chromatic polynomial of the 4-cycle graph is P(x) = 1x^4 - 4x^3 + 6x^2 - 3x.")
    print("The roots of this polynomial are:")
    
    # Print each root found
    for root in roots:
        # Use a small tolerance for checking if the imaginary part is zero
        if abs(root.imag) < 1e-9:
            print(f"- {root.real:.4f} (Real)")
        else:
            print(f"- {root:.4f} (Complex)")

    print("\nAs shown above, two of the chromatic roots are complex, which proves statement B is true.")

solve_and_check_roots()