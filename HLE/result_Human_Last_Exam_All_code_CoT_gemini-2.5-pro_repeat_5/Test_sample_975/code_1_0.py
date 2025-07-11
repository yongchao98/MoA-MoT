def print_answer_D():
    """
    This function prints the components of the magnetic field H as given in option D.
    The problem asks for the mathematical expressions for the field. We will print them
    in a structured way for clarity.
    """
    
    # Define variables for printing - these are symbolic
    M0 = "M_0"
    Rp = "R_p"
    R = "R"
    r = "r"
    cos_theta = "cos(theta)"
    sin_theta = "sin(theta)"
    i_r = "i_r"
    i_theta = "i_theta"
    
    print("Based on the derivation, option D is the correct choice, assuming a sign typo in the H-field for the inner region.")
    print("The expressions from option D are:")
    print("-" * 30)

    # Region 1: 0 < r < R_p
    H1_expression = f"H = -{M0} * (2*{Rp}^3 + {R}^3) / (3*{R}^3) * (-{cos_theta} * {i_r} + {sin_theta} * {i_theta})"
    print(f"In the region 0 < r < {Rp}:")
    print(H1_expression)
    
    # We can break down H1 into its components for clarity
    H1_r_coeff = f"-{M0} * (2*{Rp}^3 + {R}^3) / (3*{R}^3) * (-{cos_theta})"
    H1_theta_coeff = f"-{M0} * (2*{Rp}^3 + {R}^3) / (3*{R}^3) * ({sin_theta})"
    
    print("\nBreaking it down:")
    print(f"  H_r = {H1_r_coeff}")
    print(f"  H_theta = {H1_theta_coeff}")

    print("-" * 30)

    # Region 2: R_p < r < R
    H2_r_expression = f"- (2*{M0}/3) * [({Rp}/{R})^3 - ({Rp}/{r})^3] * {cos_theta} * {i_r}"
    H2_theta_expression = f"({M0}/3) * [2*({Rp}/{R})^3 + ({Rp}/{r})^3] * {sin_theta} * {i_theta}"
    
    print(f"In the region {Rp} < r < {R}:")
    print(f"H = {H2_r_expression} + {H2_theta_expression}")
    
    # We can break down H2 into its components for clarity
    H2_r_coeff = f"- (2*{M0}/3) * [({Rp}/{R})^3 - ({Rp}/{r})^3] * {cos_theta}"
    H2_theta_coeff = f"({M0}/3) * [2*({Rp}/{R})^3 + ({Rp}/{r})^3] * {sin_theta}"

    print("\nBreaking it down:")
    print(f"  H_r = {H2_r_coeff}")
    print(f"  H_theta = {H2_theta_coeff}")

print_answer_D()