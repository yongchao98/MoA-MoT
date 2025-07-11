def solve_work_equation():
    """
    This function prints the final derived formula for the work done by the current source.
    The derivation steps are explained in the text above.
    """
    
    # Define the symbols as strings for printing the final equation
    mu = "μ"
    mu_0 = "μ_0"
    N_sq = "N^2"
    w = "w"
    x2_minus_x1 = "(x_2 - x_1)"
    I2_sq_minus_I1_sq = "(I_2^2 - I_1^2)"
    two_g = "(2g)"
    
    # Print the final equation with each component clearly shown
    print("The work done by the current source for each cycle is W.")
    print("The final equation is derived as:")
    print(f"W = - (({mu} - {mu_0}) / {two_g}) * {N_sq} * {w} * {x2_minus_x1} * {I2_sq_minus_I1_sq}")

solve_work_equation()