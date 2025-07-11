import sympy
from pyknotid.make import make_knot
from pyknotid.catalogue import get_knot

def solve_knot_polynomial_difference():
    """
    This function calculates the difference in the z^2 coefficients
    of the Alexander-Conway polynomials for two knots:
    1. The closure of a given braid beta.
    2. The knot 10_4 from the Rolfsen table.
    """
    # Define the symbol for the polynomial variable
    z = sympy.Symbol('z')

    # --- Part 1: First knot from the braid beta ---

    # The braid is beta = s_4^-1 s_4^-1 s_3^-1 s_4 s_3^-1 s_2 s_1^-1 s_3^-1 s_2^-1 s_2^-1 s_2^-1 s_1^-1.
    # In pyknotid, the generator s_i is represented by the integer i,
    # and its inverse s_i^-1 is represented by -i.
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]

    # Create the knot k_beta from the closure of the braid.
    k_beta = make_knot(beta_word, braid=True)

    # Calculate the Alexander-Conway polynomial for the closure of beta.
    conway_beta = k_beta.conway_polynomial()

    # Extract the z^2 coefficient from the polynomial.
    # The .coeff(variable, power) method gets the coefficient.
    coeff_beta = conway_beta.coeff(z, 2)

    # --- Part 2: Second knot 10_4 ---

    # Get the knot 10_4 from the pyknotid catalogue.
    k_10_4 = get_knot('10_4')

    # Calculate the Alexander-Conway polynomial for 10_4.
    conway_10_4 = k_10_4.conway_polynomial()

    # Extract the z^2 coefficient from the polynomial.
    coeff_10_4 = conway_10_4.coeff(z, 2)

    # --- Part 3: Calculate the difference and present the results ---

    # Calculate the difference between the coefficients.
    difference = coeff_beta - coeff_10_4

    # Print the results in a clear format.
    # We use sympy.pretty to format the polynomial expressions nicely.
    print("The Alexander-Conway polynomial for the closure of beta is:")
    print(f"∇_β̄(z) = {sympy.pretty(conway_beta, use_unicode=False)}")
    print("\nThe Alexander-Conway polynomial for the knot 10_4 is:")
    print(f"∇_10_4(z) = {sympy.pretty(conway_10_4, use_unicode=False)}")
    print("\n--------------------------------------------------")
    print(f"The z^2 coefficient for ∇_β̄ is: {coeff_beta}")
    print(f"The z^2 coefficient for ∇_10_4 is: {coeff_10_4}")
    print("--------------------------------------------------\n")
    print("The difference in the z^2 coefficients is:")
    print(f"{coeff_beta} - {coeff_10_4} = {difference}")

if __name__ == "__main__":
    solve_knot_polynomial_difference()