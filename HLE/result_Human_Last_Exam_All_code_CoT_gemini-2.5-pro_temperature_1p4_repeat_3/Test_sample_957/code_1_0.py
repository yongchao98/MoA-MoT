def print_final_equations():
    """
    This function prints the final expressions for the overall transmission and
    reflection coefficients for the Fabry-Pérot cavity as given in the correct answer choice.
    The variables are symbolic.
    - tau: Overall transmission coefficient
    - rho: Overall reflection coefficient
    - tau_m: Transmission coefficient of a single mirror
    - rho_m: Reflection coefficient of a single mirror
    - k0: Wavenumber of the light
    - d: Distance between the mirrors
    - i: Imaginary unit
    """

    # Expressions from Option D
    tau_expression = "τ = (τ_m**2 * e**(i*k0*d)) / (1 - ρ_m**2 * e**(i*2*k0*d))"
    rho_expression = "ρ = (1 - (ρ_m - τ_m**2) * e**(i*2*k0*d) * ρ_m) / (1 - ρ_m**2 * e**(i*2*k0*d))"

    print("The final equations for the overall transmission (τ) and reflection (ρ) coefficients are:")
    print("")
    print("Transmission coefficient:")
    print(tau_expression)
    print("")
    print("Reflection coefficient:")
    print(rho_expression)

print_final_equations()