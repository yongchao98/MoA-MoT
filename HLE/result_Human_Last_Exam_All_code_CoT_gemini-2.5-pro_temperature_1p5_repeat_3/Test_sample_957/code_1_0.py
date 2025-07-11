def print_final_equations():
    """
    This function prints the final equations from the chosen answer choice (D)
    and lists the numerical constants (coefficients and exponents) as requested.
    """
    
    tau_equation = "τ = (τ_m**2 * e**(i*k_0*d)) / (1 - ρ_m**2 * e**(i*2*k_0*d))"
    rho_equation = "ρ = (1 - (ρ_m - τ_m**2) * e**(i*2*k_0*d) * ρ_m) / (1 - ρ_m**2 * e**(i*2*k_0*d))"
    
    print("The equations from the selected answer choice are:")
    print(tau_equation)
    print(rho_equation)
    
    print("\nAnalyzing the numbers in the transmission equation (τ):")
    print("Exponent on τ_m: 2")
    print("Coefficient on (i*k_0*d) in the numerator exponential: 1")
    print("Constant term in the denominator: 1")
    print("Coefficient of the second term in the denominator: -1")
    print("Exponent on ρ_m in the denominator: 2")
    print("Coefficient on (i*k_0*d) in the denominator exponential: 2")
    
    print("\nAnalyzing the numbers in the reflection equation (ρ):")
    print("Constant term in the numerator: 1")
    print("Coefficient of the second term in the numerator: -1")
    print("Exponent on ρ_m inside the parenthesis: 1")
    print("Exponent on τ_m inside the parenthesis: 2")
    print("Coefficient on (i*k_0*d) in the numerator exponential: 2")
    print("Constant term in the denominator: 1")
    print("Coefficient of the second term in the denominator: -1")
    print("Exponent on ρ_m in the denominator: 2")
    print("Coefficient on (i*k_0*d) in the denominator exponential: 2")

print_final_equations()