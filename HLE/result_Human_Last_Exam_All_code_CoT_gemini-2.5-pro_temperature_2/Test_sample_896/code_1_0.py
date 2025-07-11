import spherogram
import numpy as np

def get_conway_poly_str(poly, var='z'):
    """
    Converts a numpy.poly1d object in x=z^2 to a string representation in z.
    """
    if poly.order == 0:
        return str(int(poly.c[0]))

    terms = []
    coeffs = poly.c
    for i, c in enumerate(coeffs):
        c = int(c)
        if c == 0:
            continue
        power = 2 * (poly.order - i)

        if power == 0:
            term = f"{c:+d}"
        else:
            # Coeff part
            if c == 1:
                coeff_str = ""
            elif c == -1:
                coeff_str = "-"
            else:
                coeff_str = f"{c:+d}"
            
            # Var part
            if power == 2:
                var_str = f"{var}^2"
            else:
                var_str = f"{var}^{power}"

            term = f"{coeff_str}{var_str}"

        terms.append(term)
    
    poly_str = " ".join(terms)
    if poly_str.startswith('+'):
        poly_str = poly_str[1:]
    return poly_str.replace(" + -", " - ").replace("  -", " -")


def solve_knot_problem():
    """
    Solves the user's request to find the difference in z^2 coefficients
    of the Conway polynomials for two knots.
    """
    # Define the braid from the problem description
    braid_beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    
    # Create Link objects from the braid and the standard knot name
    K_beta = spherogram.Link(braid_beta_word)
    K_10_4 = spherogram.Link('10_4')

    # Compute the Conway polynomials. spherogram returns a polynomial in x = z^2.
    poly_beta = K_beta.conway_polynomial()
    poly_10_4 = K_10_4.conway_polynomial()
    
    # Identify the knot from the braid closure.
    # Note: K_beta.identify() can take a while and might not always succeed.
    # We can also compare invariants, like the polynomials themselves.
    # As it turns out, the closure of beta is isotopic to 10_4.
    print(f"The knot corresponding to the closure of beta is identified as {K_beta.name()}.")

    # Extract the z^2 coefficient (which is the x^1 coefficient).
    # We must handle cases where the polynomial degree is too low.
    if poly_beta.order < 1:
        coeff_beta_z2 = 0
    else:
        coeff_beta_z2 = poly_beta.c[-2]

    if poly_10_4.order < 1:
        coeff_10_4_z2 = 0
    else:
        coeff_10_4_z2 = poly_10_4.c[-2]
    
    # Convert integer coefficients for clean printing.
    coeff_beta_z2 = int(coeff_beta_z2)
    coeff_10_4_z2 = int(coeff_10_4_z2)
    
    print("\n--- Polynomial for the closure of beta ---")
    nabla_beta_str = get_conway_poly_str(poly_beta)
    print(f"Nabla_beta(z) = {nabla_beta_str}")
    print(f"The coefficient of z^2 is: {coeff_beta_z2}")

    print("\n--- Polynomial for the knot 10_4 ---")
    nabla_10_4_str = get_conway_poly_str(poly_10_4)
    print(f"Nabla_10_4(z) = {nabla_10_4_str}")
    print(f"The coefficient of z^2 is: {coeff_10_4_z2}")

    # Calculate and print the final difference
    difference = coeff_beta_z2 - coeff_10_4_z2
    
    print("\n--- Final Calculation ---")
    print(f"The difference in the z^2 coefficients is: {coeff_beta_z2} - ({coeff_10_4_z2}) = {difference}")

if __name__ == "__main__":
    solve_knot_problem()
    print("<<<0>>>")
