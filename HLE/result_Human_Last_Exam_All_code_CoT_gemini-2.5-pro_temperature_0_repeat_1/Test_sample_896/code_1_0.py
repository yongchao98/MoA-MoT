import spherogram
import sympy

def solve_knot_problem():
    """
    This function calculates the difference in the z^2 coefficients of the
    Conway polynomials for the closure of a given braid and the knot 10_4.
    """
    # Define the symbolic variable for our polynomials
    z = sympy.Symbol('z')

    # --- Step 1: Analyze the knot from the braid beta ---

    # The braid beta is in B_5, with the word:
    # s_4^-1 s_4^-1 s_3^-1 s_4 s_3^-1 s_2 s_1^-1 s_3^-1 s_2^-1 s_2^-1 s_2^-1 s_1^-1
    # We represent this numerically for spherogram.
    braid_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    
    # Create the Braid object on 5 strands
    B = spherogram.Braid(5, braid_word)

    # Create the knot K_beta by taking the closure of the braid
    K_beta = spherogram.Link(B)

    # Calculate the Conway polynomial for K_beta
    conway_poly_beta = K_beta.alexander_polynomial(variable=z, conway=True)

    # Extract the coefficient of the z^2 term
    coeff_beta = conway_poly_beta.coeff(z, 2)

    # --- Step 2: Analyze the knot 10_4 ---

    # Create the Link object for the knot 10_4
    K_10_4 = spherogram.Link('10_4')

    # Calculate the Conway polynomial for 10_4
    conway_poly_10_4 = K_10_4.alexander_polynomial(variable=z, conway=True)

    # Extract the coefficient of the z^2 term
    coeff_10_4 = conway_poly_10_4.coeff(z, 2)

    # --- Step 3: Calculate the difference and print the results ---

    difference = coeff_beta - coeff_10_4

    print(f"The Conway polynomial for the closure of beta, Nabla_beta_bar(z), is: {sympy.expand(conway_poly_beta)}")
    print(f"The Conway polynomial for 10_4, Nabla_10_4(z), is: {sympy.expand(conway_poly_10_4)}")
    print("-" * 20)
    print(f"The z^2 coefficient of Nabla_beta_bar(z) is: {coeff_beta}")
    print(f"The z^2 coefficient of Nabla_10_4(z) is: {coeff_10_4}")
    print("-" * 20)
    print(f"The difference between the coefficients is: {coeff_beta} - ({coeff_10_4}) = {difference}")

if __name__ == '__main__':
    solve_knot_problem()