import numpy as np

def solve_glyoxal_huckel():
    """
    Calculates the 4 pi-electron energy levels of glyoxal using modified Huckel theory.
    """
    # The Hückel secular determinant for glyoxal with parameters hO=1 and kO=0.8
    # leads to the characteristic polynomial:
    # x^4 + 2x^3 - 1.28x^2 - 3.28x - 0.5904 = 0
    # where E = α - xβ.

    # Coefficients of the polynomial, from the highest power to the constant term.
    coefficients = [1, 2, -1.28, -3.28, -0.5904]

    # Find the roots of the polynomial using NumPy. These are the 'x' values.
    roots_x = np.roots(coefficients)

    # The energies are E = α - xβ.
    # To list the energies from lowest to highest, we need to sort the x values
    # from largest to smallest.
    sorted_roots_x = np.sort(roots_x)[::-1]

    print("The four Huckel energies (E) for glyoxal are derived from the roots (x) of the secular polynomial.")
    print("The relationship is E = α - xβ.")
    print("\nThe calculated energies, from lowest to highest, are:")

    # Print the energy for each root in the final equation form
    for i, x_val in enumerate(sorted_roots_x):
        # Format the output string based on the sign of the root x
        if x_val < 0:
            # Equation becomes E = α + |x|β
            print(f"E{i + 1} = α + {-x_val:.4f}β")
        else:
            # Equation is E = α - xβ
            print(f"E{i + 1} = α - {x_val:.4f}β")

# Execute the function
solve_glyoxal_huckel()