import pyknotid
from pyknotid.representations import Braid
from pyknotid.knotinfo import knot_from_name
import sympy

def solve_knot_problem():
    """
    This function solves the user's request by computing and comparing
    the z^2 coefficients of the Alexander-Conway polynomials for two knots.
    """
    try:
        # Step 1: Define the braid beta from B_5
        # beta = sigma_4^{-1}sigma_4^{-1}sigma_3^{-1}sigma_4 ...
        # In pyknotid, sigma_i is represented by i and sigma_i^{-1} by -i.
        beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]

        # Step 2: Create the Braid object and compute its Conway polynomial
        print("Analyzing the knot from the closure of braid beta...")
        b = Braid(beta_word)
        k_beta = b.get_knot()
        conway_beta = k_beta.conway_polynomial

        # Step 3: Extract the z^2 coefficient for beta_bar
        # The polynomial variable (usually z) is found from the expression
        poly_vars_beta = conway_beta.free_symbols
        if not poly_vars_beta:
             # Handle constant polynomials (like for the unknot)
            z_var_beta = sympy.Symbol('z')
        else:
            z_var_beta = list(poly_vars_beta)[0]
        
        coeff_beta = conway_beta.coeff(z_var_beta, 2)
        print(f"The Alexander-Conway polynomial for the closure of beta is: {sympy.pretty(conway_beta)}")
        print(f"The z^2 coefficient is: {coeff_beta}")
        print("-" * 20)

        # Step 4: Get the knot 10_4 and its Conway polynomial
        print("Analyzing the knot 10_4...")
        k_10_4 = knot_from_name('10_4')
        conway_10_4 = k_10_4.conway_polynomial
        
        # Step 5: Extract the z^2 coefficient for 10_4
        poly_vars_10_4 = conway_10_4.free_symbols
        z_var_10_4 = list(poly_vars_10_4)[0]
        coeff_10_4 = conway_10_4.coeff(z_var_10_4, 2)
        print(f"The Alexander-Conway polynomial for the knot 10_4 is: {sympy.pretty(conway_10_4)}")
        print(f"The z^2 coefficient is: {coeff_10_4}")
        print("-" * 20)

        # Step 6: Calculate and display the difference
        difference = coeff_beta - coeff_10_4
        print("The difference in the z^2 coefficients is:")
        print(f"(coefficient for beta_bar) - (coefficient for 10_4) = ({coeff_beta}) - ({coeff_10_4}) = {difference}")

    except ImportError:
        print("This script requires the 'pyknotid' and 'sympy' libraries.")
        print("Please install them by running: pip install pyknotid sympy")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    solve_knot_problem()