def print_fabry_perot_coefficients():
    """
    Prints the expressions for the overall reflection and transmission coefficients
    for a Fabry-Pérot cavity based on the derived result.
    """
    
    # The expressions correspond to option D.
    # tau is the overall transmission coefficient.
    # rho is the overall reflection coefficient.
    # tau_m is the transmission coefficient of a single mirror.
    # rho_m is the reflection coefficient of a single mirror.
    # k0 is the wave number in free space.
    # d is the distance between the mirrors.
    
    tau_expression = "tau = (tau_m**2 * exp(i*k0*d)) / (1 - rho_m**2 * exp(i*2*k0*d))"
    rho_expression = "rho = (1 - (rho_m - tau_m**2) * exp(i*2*k0*d) * rho_m) / (1 - rho_m**2 * exp(i*2*k0*d))"
    
    print("The overall transmission coefficient (tau) is:")
    print(tau_expression)
    print("\nThe overall reflection coefficient (rho) is:")
    print(rho_expression)

print_fabry_perot_coefficients()