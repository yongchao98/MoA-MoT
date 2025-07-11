def print_fabry_perot_formulas():
    """
    This function prints the derived formulas for the overall transmission (tau)
    and reflection (rho) coefficients of a Fabry-Pérot cavity.
    """

    # Define the variables symbolically as strings
    tau_m = "τ_m"
    rho_m = "ρ_m"
    k0 = "k₀"
    d = "d"
    i = "i" # Imaginary unit

    # Construct the formulas as strings
    # Note: In Python code, complex numbers are usually written with 'j', but 'i' is common in physics.
    # We use ** for exponentiation and * for multiplication for clarity.
    
    # Transmission coefficient formula
    tau_numerator = f"{tau_m}**2 * exp({i}*{k0}*{d})"
    tau_denominator = f"1 - {rho_m}**2 * exp({i}*2*{k0}*{d})"
    tau_formula = f"({tau_numerator}) / ({tau_denominator})"

    # Reflection coefficient formula
    rho_numerator_part2 = f"{tau_m}**2 * {rho_m} * exp({i}*2*{k0}*{d})"
    rho_denominator = f"1 - {rho_m}**2 * exp({i}*2*{k0}*{d})"
    rho_formula = f"{rho_m} + ({rho_numerator_part2}) / ({rho_denominator})"

    # Print the results
    print("Based on the summation of multiple beams, the derived coefficients are:")
    print("-" * 30)
    
    print("Overall Transmission Coefficient (τ):")
    print(f"τ = {tau_formula}")
    print("\nThis expression for τ matches option D.")
    
    print("\nOverall Reflection Coefficient (ρ):")
    print(f"ρ = {rho_formula}")
    print("\nThis is the correct expression for ρ under the simplified model.")
    print("The expression for τ in option D is correct, so D is the intended answer, despite an apparent typo in its expression for ρ.")

# Execute the function to print the formulas
print_fabry_perot_formulas()