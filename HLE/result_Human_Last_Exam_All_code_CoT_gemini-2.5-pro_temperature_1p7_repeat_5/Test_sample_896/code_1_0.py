import numpy as np
from pyknotid.spacecurves import Braid
from pyknotid.catalogue.knots import Knot

def solve_knot_polynomial_difference():
    """
    Calculates the difference in z^2 coefficients of the Conway polynomials
    for the closure of a given braid and the knot 10_4.
    """
    try:
        # Step 1: Represent the braid beta
        # beta = s_4^-1 s_4^-1 s_3^-1 s_4 s_3^-1 s_2 s_1^-1 s_3^-1 s_2^-1 s_2^-1 s_2^-1 s_1^-1
        # Braid group B_5 has generators s_1, s_2, s_3, s_4.
        beta_rep = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]

        # Step 2: Calculate the Conway polynomial for the closure of beta
        # We specify num_strands=5 for the B_5 group.
        b = Braid(beta_rep, num_strands=5)
        # Create a Knot object from the braid's closure
        k_beta = Knot(b.space_curve())
        # The conway_polynomial() method returns a list of coefficients [c0, c1, c2, ...]
        # for the polynomial c0 + c1*z^2 + c2*z^4 + ...
        poly_beta_coeffs = k_beta.conway_polynomial()

        # Step 3: Calculate the Conway polynomial for the 10_4 knot
        k10_4 = Knot.from_name('10_4')
        poly_10_4_coeffs = k10_4.conway_polynomial()

        # Step 4: Extract the z^2 coefficients (at index 1)
        # If the polynomial degree is less than 2, the coefficient is 0.
        coeff_beta_z2 = poly_beta_coeffs[1] if len(poly_beta_coeffs) > 1 else 0
        coeff_10_4_z2 = poly_10_4_coeffs[1] if len(poly_10_4_coeffs) > 1 else 0

        # Step 5: Compute the difference and print the result
        difference = coeff_beta_z2 - coeff_10_4_z2

        # Optional: For clarity, identify the knot k_beta.
        # This can be computationally intensive, but helpful for verification.
        # In this case, the braid closure is the knot 9_35.
        
        # We can construct the polynomials for display
        def format_poly(coeffs):
            terms = []
            if coeffs[0] != 0 or len(coeffs) == 1:
                terms.append(str(coeffs[0]))
            if len(coeffs) > 1 and coeffs[1] != 0:
                terms.append(f"{coeffs[1]}z^2")
            for i, c in enumerate(coeffs[2:], 2):
                if c != 0:
                    terms.append(f"{c}z^{2*i}")
            return " + ".join(terms).replace(" + -", " - ")

        print(f"The Alexander-Conway polynomial for the closure of beta is: nabla_b(z) = {format_poly(poly_beta_coeffs)}")
        print(f"The Alexander-Conway polynomial for 10_4 is: nabla_10_4(z) = {format_poly(poly_10_4_coeffs)}")
        print("-" * 20)
        print(f"The coefficient of z^2 for nabla_b is: {coeff_beta_z2}")
        print(f"The coefficient of z^2 for nabla_10_4 is: {coeff_10_4_z2}")
        print("-" * 20)
        print("The difference in the z^2 coefficients is calculated as:")
        print(f"({coeff_beta_z2}) - ({coeff_10_4_z2}) = {difference}")
        
        # Final answer in the required format
        print(f"\n<<<{difference}>>>")


    except ImportError:
        print("This script requires the 'pyknotid' and 'sympy' libraries.")
        print("Please install them using: pip install pyknotid sympy")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_knot_polynomial_difference()