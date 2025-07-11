import numpy as np

def solve_coordinates_average():
    """
    This function finds the complex coordinates z where the inverse of the gradient of S(z)
    approaches zero and calculates their average.

    The coordinates are the poles of the B(z) field, determined by the roots of the
    denominator in its governing equation.
    """
    # The relevant polynomial is 4*z^4 - z^3 + z^2 + 1 = 0.
    # The coefficients are for z^4, z^3, z^2, z^1, z^0 respectively.
    coeffs = [4, -1, 1, 0, 1]

    # Find the roots of the polynomial.
    roots = np.roots(coeffs)

    # Calculate the average of the roots.
    average_of_roots = np.mean(roots)

    print("The complex coordinates z are the roots of the equation 4z^4 - z^3 + z^2 + 1 = 0.")
    print("The identified coordinates are:")
    for i, root in enumerate(roots):
        print(f"z{i+1} = {root}")

    # Display the final calculation as an equation
    print("\nThe average of these coordinates is calculated as:")
    equation_str = f"({roots[0]}) + ({roots[1]}) + ({roots[2]}) + ({roots[3]}) / 4"
    print(f"{equation_str} = {average_of_roots}")
    
    # In this case, the imaginary parts cancel out and we are left with a real number.
    print(f"\nThe sum of roots from Vieta's formulas is -(-1)/4 = 0.25.")
    print(f"The average of roots is (0.25) / 4 = {0.25/4}.")


solve_coordinates_average()