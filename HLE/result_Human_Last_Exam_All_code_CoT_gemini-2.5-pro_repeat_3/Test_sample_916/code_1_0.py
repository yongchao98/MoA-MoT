def print_force_equation():
    """
    Prints the symbolic formula for the instantaneous force f_x(t).
    """
    # Define symbolic variables for clarity
    term_constant = "-2*pi*R*N"
    term_numerator = "mu_0 * (1 - alpha_T * (T - T_0)) * N_0 * I_0 * i_0 * sin(omega*t)"
    term_denominator_g = "g^2"
    term_denominator_saturation = "(1 + (mu_0 * N_0 * I_0)/(g * B_s))"
    
    # Assemble and print the equation for f_x(t)
    print("f_x(t) = {} * (({}) / ({} * {}))".format(
        term_constant, 
        term_numerator, 
        term_denominator_g, 
        term_denominator_saturation
    ))
    
    # Print the equation with each term identified as per the prompt's spirit.
    print("\nWhere the equation is:")
    print("f_x(t) = -2 * pi * R * N * (mu_0 * (1 - alpha_T * (T - T_0)) * N_0 * I_0 * i_0 * sin(omega*t)) / (g^2 * (1 + (mu_0 * N_0 * I_0)/(g * B_s)))")

print_force_equation()