def print_fabry_perot_coefficients():
    """
    Prints the expressions for the overall reflection and transmission 
    coefficients of the Fabry-Pérot cavity based on the derived answer.
    """
    
    # The derivation of the transmission coefficient τ by summing the geometric series
    # of multiply-reflected waves yields a result that matches Option D.
    # While the expression for ρ in the options is non-standard, the correctness
    # of τ makes Option D the correct choice.
    
    # Define the variables for clarity in the final equations
    tau_m_squared = "τ_m^2"
    rho_m_squared = "ρ_m^2"
    rho_m = "ρ_m"
    phase_factor_d = "e^(i*k_0*d)"
    phase_factor_2d = "e^(i*2*k_0*d)"

    # Transmission coefficient (τ) from Option D
    tau_numerator = f"{tau_m_squared} * {phase_factor_d}"
    tau_denominator = f"1 - {rho_m_squared} * {phase_factor_2d}"
    
    # Reflection coefficient (ρ) from Option D
    rho_numerator = f"1 - ({rho_m} - {tau_m_squared}) * {phase_factor_2d} * {rho_m}"
    rho_denominator = f"1 - {rho_m_squared} * {phase_factor_2d}"

    print("The overall transmission coefficient τ is given by the equation:")
    print(f"τ = ({tau_numerator}) / ({tau_denominator})")
    
    print("\nThe overall reflection coefficient ρ is given by the equation:")
    print(f"ρ = ({rho_numerator}) / ({rho_denominator})")

# Execute the function to display the results
print_fabry_perot_coefficients()