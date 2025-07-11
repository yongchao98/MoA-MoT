import sympy
from pyknotid.representations import Braid
from pyknotid.catalogue import get_knot

def solve_knot_polynomial_coefficient_difference():
    """
    This script calculates the difference in the z^2 coefficients of the
    Alexander-Conway polynomials for two knots: one given as a braid closure
    and the other as the knot 10_4.
    
    It requires the 'pyknotid' library. You can install it via pip:
    pip install pyknotid
    """
    
    # Define the symbolic variable for the polynomial
    z = sympy.Symbol('z')

    # 1. Define the braid beta from the problem description.
    # The braid word is given in terms of generators s_i.
    # In pyknotid, s_i is represented by the integer i, and its inverse s_i^-1 by -i.
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    num_strands = 5
    
    try:
        # Create the Braid object
        braid_beta = Braid(braid_word=beta_word, num_strands=num_strands)

        # 2. Calculate the Conway polynomial for the closure of beta
        nabla_beta = braid_beta.conway_polynomial
        
        # 3. Extract the coefficient of z^2 from the polynomial
        # Convert to a sympy.Poly object to safely extract coefficients
        nabla_beta_poly = sympy.Poly(nabla_beta, z)
        coeff_beta = nabla_beta_poly.coeff_monomial(z**2)

        # 4. Get the knot 10_4 from the pyknotid catalogue
        knot_10_4 = get_knot('10_4')
        
        # 5. Get its Conway polynomial
        nabla_10_4 = knot_10_4.conway_polynomial
        
        # 6. Extract the coefficient of z^2
        nabla_10_4_poly = sympy.Poly(nabla_10_4, z)
        coeff_10_4 = nabla_10_4_poly.coeff_monomial(z**2)
        
        # 7. Compute the difference and print the final equation
        difference = coeff_beta - coeff_10_4
        
        print(f"The Alexander-Conway polynomial for the closure of the braid beta is: {nabla_beta}")
        print(f"The Alexander-Conway polynomial for the knot 10_4 is: {nabla_10_4}")
        print(f"The coefficient of z^2 in nabla_beta is: {coeff_beta}")
        print(f"The coefficient of z^2 in nabla_10_4 is: {coeff_10_4}")
        print(f"The difference between the coefficients is: {coeff_beta} - {coeff_10_4} = {difference}")

    except ImportError:
        print("Please install the required libraries 'pyknotid', 'sympy', and 'numpy'.")
        print("You can install them using: pip install pyknotid")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_knot_polynomial_coefficient_difference()