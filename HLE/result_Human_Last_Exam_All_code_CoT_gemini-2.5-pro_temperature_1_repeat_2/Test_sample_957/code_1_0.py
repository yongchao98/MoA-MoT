def print_fabry_perot_coefficients():
    """
    This function prints the final equations for the overall reflection and transmission
    coefficients of the Fabry-Pérot cavity based on the chosen answer choice.
    """
    # Variables are represented as strings for display purposes.
    tau_m = "τ_m"
    rho_m = "ρ_m"
    k0d = "k₀d"

    # The formulas from the selected answer choice D
    
    # Transmission coefficient τ
    tau_numerator = f"{tau_m}² * exp(i*{k0d})"
    tau_denominator = f"1 - {rho_m}² * exp(i*2*{k0d})"
    
    # Reflection coefficient ρ
    rho_numerator = f"1 - ({rho_m} - {tau_m}²) * exp(i*2*{k0d}) * {rho_m}"
    rho_denominator = f"1 - {rho_m}² * exp(i*2*{k0d})"

    print("The overall transmission and reflection coefficients are:")
    print("-" * 50)
    print(f"Transmission coefficient (τ):")
    print(f"      {tau_numerator}")
    print(f"τ = -----------------------------")
    print(f"      {tau_denominator}")
    print("\n")
    print(f"Reflection coefficient (ρ):")
    print(f"      {rho_numerator}")
    print(f"ρ = -----------------------------")
    print(f"      {rho_denominator}")
    print("-" * 50)

# Execute the function to display the result
print_fabry_perot_coefficients()