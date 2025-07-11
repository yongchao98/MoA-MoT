def polymer_force_law():
    """
    This function prints the derived force law for a thermally isolated
    freely jointed polymer chain.
    """
    # Define string representations of the variables for clarity
    var_x = "x"
    var_ell = "l"
    var_n = "n"
    var_E0 = "E(0)"
    
    # The numerical constants present in the final equation
    number_3_coeff = 3
    number_2_pow_n = 2
    number_2_pow_ell = 2
    
    number_3_exp = 3
    number_2_exp_x = 2
    number_2_exp_den = 2

    # Construct the final equation string piece by piece
    coefficient = f"({number_3_coeff} * {var_E0}) / ({var_n}**{number_2_pow_n} * {var_ell}**{number_2_pow_ell})"
    exponent_argument = f"-({number_3_exp} * {var_x}**{number_2_exp_x}) / ({number_2_exp_den} * {var_n}**{number_2_pow_n} * {var_ell}**{number_2_pow_ell})"
    force_equation = f"f({var_x}) = - {coefficient} * {var_x} * exp({exponent_argument})"

    # Print the final result
    print("The force law between the polymer ends under thermal isolation is:")
    print(force_equation)
    print("\nWhere:")
    print("  f(x) is the attractive force between the ends.")
    print("  x is the separation distance between the ends.")
    print("  n is the number of links in the polymer chain (assumed large).")
    print("  l is the length of a single link.")
    print("  E(0) is the total kinetic energy of the polymer at zero extension (x=0).")

polymer_force_law()
<<<f(x) = - (3 * E(0) / (n**2 * l**2)) * x * exp(-(3 * x**2) / (2 * n**2 * l**2))>>>