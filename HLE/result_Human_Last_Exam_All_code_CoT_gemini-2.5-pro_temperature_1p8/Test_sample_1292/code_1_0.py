def print_field_equation():
    """
    Prints the terms of the derived field equation for Symmetric Teleparallel Gravity.
    
    The equation represents the variation of the given action with respect to the metric.
    Each term is printed on a new line for clarity.
    """
    
    # Constants
    kappa = "8 * pi * G / c**4"
    
    # Define the terms of the equation as strings
    term1_str = "-2/sqrt(-g) * d_alpha(sqrt(-g) * P^alpha_munu)"
    term2_str = "- 2 * P_mu,alpha,beta * Q_nu^alpha,beta"
    term3_str = "+ Q^alpha,beta_mu * P_alpha,beta,nu"
    term4_str = "- 1/2 * Q * g_munu"
    rhs_str = f"= ({kappa}) * T_munu"
    
    # Print the equation
    print("The derived field equation is:")
    print(f"{term1_str} \\")
    print(f"    {term2_str} \\")
    print(f"    {term3_str} \\")
    print(f"    {term4_str} \\")
    print(f"    {rhs_str}")
    
    print("\nSymbolically, this corresponds to:")
    equation_string = (
        f"  -2/sqrt(-g) * partial_alpha(sqrt(-g)*P^alpha_munu) "
        f"- 2*P_mu_alpha_beta*Q_nu^alpha_beta "
        f"+ Q^alpha_beta_mu*P_alpha_beta_nu "
        f"- (1/2)*Q*g_munu = ({kappa})*T_munu"
    )
    print(equation_string)

    print("\nRepresenting each term from Choice B with its numerical coefficient:")
    # Mapping to choice B
    print("- coefficient of first term: -2")
    print("- coefficient of second term: -2")
    print("- coefficient of third term: 1") # +1
    print("- coefficient of fourth term: -1/2")
    
print_field_equation()