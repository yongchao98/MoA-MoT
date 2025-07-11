def print_fabry_perot_coefficients():
    """
    Prints the expressions for the overall reflection and transmission coefficients
    for a Fabry-Pérot cavity based on the selected answer choice.
    """
    # These are the formulas from option D.
    transmission_coeff_str = "τ = (τ_m^2 * e^(i*k_0*d)) / (1 - ρ_m^2 * e^(i*2*k_0*d))"
    reflection_coeff_str = "ρ = (1 - (ρ_m - τ_m^2) * e^(i*2*k_0*d) * ρ_m) / (1 - ρ_m^2 * e^(i*2*k_0*d))"

    print("The overall transmission and reflection coefficients for the Fabry-Pérot cavity are:")
    print(f"Transmission coefficient: {transmission_coeff_str}")
    print(f"Reflection coefficient: {reflection_coeff_str}")

print_fabry_perot_coefficients()