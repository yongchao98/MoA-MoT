def print_fabry_perot_coefficients():
    """
    This function prints the overall reflection and transmission coefficients for a Fabry-Pérot cavity.
    The equations correspond to the correct choice based on physics derivation.
    """
    
    # Define the expressions for the transmission (tau) and reflection (rho) coefficients
    # based on the analysis of the Fabry-Pérot cavity.
    # tau_m: single mirror transmission coefficient
    # rho_m: single mirror reflection coefficient
    # k0: wave number
    # d: distance between mirrors
    
    tau_expression = "τ = (τ_m^2 * e^(i*k_0*d)) / (1 - ρ_m^2 * e^(i*2*k_0*d))"
    rho_expression = "ρ = (1 - (ρ_m - τ_m^2) * e^(i*2*k_0*d) * ρ_m) / (1 - ρ_m^2 * e^(i*2*k_0*d))"

    print("The overall transmission and reflection coefficients for the Fabry-Pérot cavity are:")
    print("\nTransmission Coefficient (τ):")
    print(tau_expression)
    
    print("\nReflection Coefficient (ρ):")
    print(rho_expression)

# Execute the function to print the results
print_fabry_perot_coefficients()