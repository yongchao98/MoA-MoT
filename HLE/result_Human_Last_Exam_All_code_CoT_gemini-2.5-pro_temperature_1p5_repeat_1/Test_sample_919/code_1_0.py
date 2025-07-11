def solve_emi_force():
    """
    This function prints the derived formula for the force per unit area
    on the conducting plane.
    """
    # Define the components of the formula as strings
    force_vector_name = "f_vec/area"
    
    # Numerator part
    numerator_factor = "- (mu_0 / 2)"
    numerator_terms = "K_0^2 * sin^2(ay)"
    
    # Denominator part
    denominator_terms = "[cosh(ad) + (mu_0 / mu) * sinh(ad)]^2"
    
    # Direction vector
    direction = "i_x_hat"

    # Assemble the final equation string
    final_equation = f"{force_vector_name} = {numerator_factor} * ({numerator_terms}) / ({denominator_terms}) * {direction}"

    print("The derived force per unit area on the x=d interface is:")
    print(final_equation)
    
    print("\nThis corresponds to the expression:")
    print("f/area = - (mu_0 / 2) * (K_0^2 * sin^2(ay)) / (cosh(ad) + (mu_0 / mu) * sinh(ad))^2 * i_x")


solve_emi_force()