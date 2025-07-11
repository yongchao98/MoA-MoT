def print_cavity_equations():
    """
    This function prints the final equations for the overall reflection and
    transmission coefficients of a Fabry-Pérot cavity based on the derived answer.
    """
    # Define string representations for the variables for clear output
    tau_m_str = "τ_m"
    rho_m_str = "ρ_m"
    k0d_str = "k₀d" # Using subscript zero for clarity

    # Construct the string for the transmission coefficient τ from option D
    tau_numerator = f"{tau_m_str}² * exp(i*{k0d_str})"
    tau_denominator = f"1 - {rho_m_str}² * exp(i*2*{k0d_str})"
    tau_eq = f"τ = ({tau_numerator}) / ({tau_denominator})"

    # Construct the string for the reflection coefficient ρ from option D
    # Note: As determined in the derivation, this expression likely contains a typo
    # compared to the standard derivation, but we present the one from the answer choice.
    rho_numerator = f"1 - ({rho_m_str} - {tau_m_str}²) * exp(i*2*{k0d_str}) * {rho_m_str}"
    rho_denominator = f"1 - {rho_m_str}² * exp(i*2*{k0d_str})"
    rho_eq = f"ρ = ({rho_numerator}) / ({rho_denominator})"

    # Print the final equations
    print("The overall transmission coefficient τ is:")
    print(tau_eq)
    print("\nThe overall reflection coefficient ρ is:")
    print(rho_eq)

# Execute the function to display the result
print_cavity_equations()