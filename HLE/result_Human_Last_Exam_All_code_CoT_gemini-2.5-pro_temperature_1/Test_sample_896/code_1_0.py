import sympy
from pyknotid.braids import Braid
from pyknotid.catalogue import get_knot

def solve_knot_problem():
    """
    This function solves the given problem by calculating and comparing
    the z^2 coefficients of the Alexander-Conway polynomials for two knots.
    """
    # The braid β is given as σ_4⁻¹σ_4⁻¹σ_3⁻¹σ_4 σ_3⁻¹σ_2σ_1⁻¹σ_3⁻¹σ_2⁻¹σ_2⁻¹σ_2⁻¹σ_1⁻¹.
    # In pyknotid, σ_i is represented by the integer i, and its inverse σ_i⁻¹ by -i.
    # The braid group B_5 has generators σ_1, σ_2, σ_3, σ_4.
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    
    # Create a Braid object for β, specifying it has 5 strands.
    braid_beta = Braid(beta_word, num_strands=5)
    
    # Get the knot K = β-bar by taking the closure of the braid.
    knot_from_beta = braid_beta.get_knot()
    
    # Calculate the Alexander-Conway polynomial for the knot from β.
    poly_beta = knot_from_beta.conway_polynomial()
    
    # Get the reference knot 10_4 from the pyknotid catalogue.
    knot_10_4 = get_knot('10_4')
    
    # Calculate the Alexander-Conway polynomial for 10_4.
    poly_10_4 = knot_10_4.conway_polynomial()
    
    # Define the symbolic variable 'z' to work with the polynomials.
    z = sympy.var('z')
    
    # Extract the coefficient of z^2 from both polynomials.
    # The .coeff(z, 2) method returns 0 if the term doesn't exist.
    coeff_beta_z2 = poly_beta.coeff(z, 2)
    coeff_10_4_z2 = poly_10_4.coeff(z, 2)
    
    # Calculate the difference between the two coefficients.
    difference = coeff_beta_z2 - coeff_10_4_z2
    
    # Print the results in a clear, step-by-step format.
    print(f"The Alexander-Conway polynomial for the closure of β is: ∇_β(z) = {poly_beta}")
    print(f"The Alexander-Conway polynomial for 10_4 is: ∇_10_4(z) = {poly_10_4}")
    print(f"The z^2 coefficient for ∇_β(z) is: {coeff_beta_z2}")
    print(f"The z^2 coefficient for ∇_10_4(z) is: {coeff_10_4_z2}")
    print(f"The difference between the z^2 coefficients is: {coeff_beta_z2} - ({coeff_10_4_z2}) = {difference}")

if __name__ == "__main__":
    solve_knot_problem()
