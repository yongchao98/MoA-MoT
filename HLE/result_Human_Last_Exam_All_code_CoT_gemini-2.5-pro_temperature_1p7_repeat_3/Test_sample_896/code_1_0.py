# This script requires the 'pyknotid' and 'sympy' libraries.
# You can install them using pip:
# pip install pyknotid sympy

import sympy
from pyknotid.catalogue import from_braid
from pyknotid.catalogue.knots import Knot

def solve_knot_problem():
    """
    Calculates the difference in the z^2 coefficients of the Conway
    polynomial for the closure of a given braid beta and the knot 10_4.
    """
    # Define the braid word for β based on the problem description:
    # β = σ_4⁻¹σ_4⁻¹σ_3⁻¹σ_4 σ_3⁻¹σ_2σ_1⁻¹σ_3⁻¹σ_2⁻¹σ_2⁻¹σ_2⁻¹σ_1⁻¹
    # In numerical notation, σ_i is i and σ_i⁻¹ is -i.
    # The braid group is B_5, so generators are 1, 2, 3, 4.
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]

    # Create a knot object from the closure of the braid β
    k_beta = from_braid(beta_word)

    # Calculate the Conway polynomial for the closure of β
    poly_beta = k_beta.conway_polynomial()

    # Create a knot object for 10_4 from the Rolfsen catalogue
    k10_4 = Knot(10, 4)
    
    # Calculate the Conway polynomial for 10_4
    poly_10_4 = k10_4.conway_polynomial()
    
    # Define the symbolic variable z
    z = sympy.Symbol('z')

    # Extract the coefficient of z^2 from both polynomials
    # The .coeff method returns the coefficient of a given term.
    # If the term does not exist, it returns 0.
    coeff_beta = sympy.sympify(poly_beta).coeff(z, 2)
    coeff_10_4 = sympy.sympify(poly_10_4).coeff(z, 2)

    # Calculate the difference between the coefficients
    difference = coeff_beta - coeff_10_4

    # Print the results as requested
    print(f"The Alexander-Conway polynomial for the closure of the braid is: {poly_beta}")
    print(f"The Alexander-Conway polynomial for the knot 10_4 is: {poly_10_4}")
    print(f"The z^2 coefficient for the braid closure is: {coeff_beta}")
    print(f"The z^2 coefficient for 10_4 is: {coeff_10_4}")
    print(f"The difference is {coeff_beta} - {coeff_10_4} = {difference}")
    
    # Returning the final numerical answer for the submission format
    return difference

if __name__ == '__main__':
    final_difference = solve_knot_problem()
    # The problem asks for the answer in a special format.
    # Although the print statements above show the work,
    # this print statement will provide the final answer for easy parsing.
    # print(f"<<<{final_difference}>>>")
    
# We will directly print the final value as the main output of the script execution.
# This structure with a main function is good practice, but for the final output, 
# let's just make sure the printouts are clear as requested.
# I will call the function directly to produce the output.

solve_knot_problem()