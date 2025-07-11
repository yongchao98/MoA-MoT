def print_fabry_perot_coefficients():
    """
    Prints the expressions for the overall reflection and transmission coefficients
    for a Fabry-Pérot cavity as given in the selected answer choice.
    """
    
    # Define the variables as strings for printing the formula
    tau_m = "τ_m"
    rho_m = "ρ_m"
    k0d = "k_0*d"
    
    # Expression for the transmission coefficient τ
    tau_expression = f"τ = ({tau_m}^2 * e^(i*{k0d})) / (1 - {rho_m}^2 * e^(i*2*{k0d}))"
    
    # Expression for the reflection coefficient ρ
    rho_expression = f"ρ = (1 - ({rho_m} - {tau_m}^2) * e^(i*2*{k0d}) * {rho_m}) / (1 - {rho_m}^2 * e^(i*2*{k0d}))"

    # Print the final equations
    print("The overall transmission and reflection coefficients are:")
    print(tau_expression)
    print(rho_expression)

print_fabry_perot_coefficients()