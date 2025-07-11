import sympy as sp

def display_fabry_perot_coefficients():
    """
    Displays the expressions for the overall reflection and transmission coefficients
    for a Fabry-Pérot cavity as given in the chosen answer choice.
    """
    # Define the symbols
    tau_m, rho_m, k0, d = sp.symbols('tau_m rho_m k_0 d')
    I = sp.I

    # Phase factor for one pass
    phi = k0 * d

    # Expressions from Option D
    tau_expression = (tau_m**2 * sp.exp(I * phi)) / (1 - rho_m**2 * sp.exp(I * 2 * phi))
    rho_expression = (1 - (rho_m - tau_m**2) * sp.exp(I * 2 * phi) * rho_m) / (1 - rho_m**2 * sp.exp(I * 2 * phi))
    
    # Print the formulas in a readable format
    print("According to the selected option:")
    print("\nOverall Transmission Coefficient (τ):")
    sp.pretty_print(tau_expression)
    
    print("\nOverall Reflection Coefficient (ρ):")
    sp.pretty_print(rho_expression)

if __name__ == '__main__':
    display_fabry_perot_coefficients()