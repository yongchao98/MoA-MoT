def print_fabry_perot_coefficients():
    """
    This function prints the final equations for the overall reflection and
    transmission coefficients of the Fabry-Pérot cavity as given in the correct answer choice.
    """
    
    # The expressions are taken from the correct answer choice.
    # Note: e^(...) represents the exponential function, i is the imaginary unit.
    
    transmission_coefficient_formula = "τ = (τ_m^2 * e^(i*k_0*d)) / (1 - ρ_m^2 * e^(i*2*k_0*d))"
    reflection_coefficient_formula = "ρ = (1 - (ρ_m - τ_m^2) * e^(i*2*k_0*d) * ρ_m) / (1 - ρ_m^2 * e^(i*2*k_0*d))"

    print("The derived overall transmission coefficient is:")
    print(transmission_coefficient_formula)
    
    print("\nThe derived overall reflection coefficient is:")
    print(reflection_coefficient_formula)

# Execute the function to print the result.
print_fabry_perot_coefficients()