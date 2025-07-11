import spherogram
import sympy

def solve_knot_polynomial_difference():
    """
    This function calculates the difference in the z^2 coefficients of the
    Alexander-Conway polynomials for a given braid closure and the knot 10_4.
    """
    # Define the braid word for beta from the problem description.
    # In spherogram, for B_n, generators are 1, 2, ..., n-1.
    # The problem uses sigma_1, ..., sigma_4 for B_5.
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]

    # Create the Braid object for B_5
    beta = spherogram.Braid(5, beta_word)

    # Get the link from the closure of the braid
    link_beta_closure = beta.link()

    # Calculate the Alexander-Conway polynomial for the closure of beta
    poly_beta_closure = link_beta_closure.alexander_conway_polynomial()

    # Define the symbolic variable z to extract coefficients
    z = sympy.Symbol('z')

    # Extract the z^2 coefficient from the polynomial of beta's closure
    coeff_beta_closure = poly_beta_closure.coeff(z**2)

    print(f"The Alexander-Conway polynomial for the closure of beta is: {poly_beta_closure}")
    print(f"The z^2 coefficient for the closure of beta is: {coeff_beta_closure}")
    print("-" * 30)

    # Create the Link object for the knot 10_4
    link_10_4 = spherogram.Link('10_4')

    # Calculate the Alexander-Conway polynomial for 10_4
    poly_10_4 = link_10_4.alexander_conway_polynomial()

    # Extract the z^2 coefficient from the polynomial of 10_4
    coeff_10_4 = poly_10_4.coeff(z**2)

    print(f"The Alexander-Conway polynomial for the knot 10_4 is: {poly_10_4}")
    print(f"The z^2 coefficient for 10_4 is: {coeff_10_4}")
    print("-" * 30)

    # Calculate the difference between the coefficients
    difference = coeff_beta_closure - coeff_10_4

    # Print the final equation as requested
    print("The difference in the z^2 coefficients is:")
    print(f"{coeff_beta_closure} - ({coeff_10_4}) = {difference}")

if __name__ == '__main__':
    solve_knot_polynomial_difference()