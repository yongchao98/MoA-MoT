def calculate_work_formula():
    """
    This function prints the derived formula for the work done by the current source.
    
    The variables in the formula are:
    μ (mu): Magnetic permeability of the movable block
    μ₀ (mu_0): Magnetic permeability of the air gap
    N: Number of turns in the winding
    w: Width of the cross-sectional area
    g: Length of the air gap
    x₁, x₂: Initial and final positions of the block
    I₁, I₂: Initial and final currents
    """

    # Representing the components of the final formula as strings
    mu = "μ"
    mu_0 = "μ₀"
    N = "N"
    w_var = "w"
    g_var = "g"
    x1 = "x₁"
    x2 = "x₂"
    I1 = "I₁"
    I2 = "I₂"
    
    # Constructing the final formula string
    # W = - (μ - μ₀) / (2*g) * N² * w * (x₂ - x₁) * (I₂² - I₁²)
    numerator_part1 = f"({mu} - {mu_0})"
    numerator_part2 = f"{N}²"
    numerator_part3 = f"{w_var}"
    numerator_part4 = f"({x2} - {x1})"
    numerator_part5 = f"({I2}² - {I1}²)"
    denominator = f"2{g_var}"

    print(f"The formula for the work (W) done by the current source per cycle is:")
    # Using '−' for the minus sign for better readability
    print(f"W = − ( {numerator_part1} * {numerator_part2} * {numerator_part3} * {numerator_part4} * {numerator_part5} ) / ( {denominator} )")
    
    # For verification, here is the formula with variable names used in code:
    mu = "mu"
    mu_0 = "mu_0"
    w = "w"
    g = "g"
    x_1 = "x_1"
    x_2 = "x_2"
    I_1 = "I_1"
    I_2 = "I_2"
    print("\nIn standard programming notation:")
    # Printing each term as requested.
    # The numbers in the final equation are 2, 2, 2, 1
    term1 = f"({mu} - {mu_0})"
    term2 = f"{N}**2"
    term3 = f"{w}"
    term4 = f"({x_2} - {x_1})"
    term5 = f"({I_2}**2 - {I_1}**2)"
    denom_term = f"2*{g}"
    print(f"W = - ( {term1} * {term2} * {term3} * {term4} * {term5} ) / ( {denom_term} )")

calculate_work_formula()