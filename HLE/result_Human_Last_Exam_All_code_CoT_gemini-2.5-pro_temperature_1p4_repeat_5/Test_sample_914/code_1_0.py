import sympy

def print_force_equation():
    """
    This function prints the derived symbolic equation for the force F_x.
    """
    # Define symbolic variables for clarity in the output.
    a, D, mu_0, I_0, sigma_1, sigma_2 = sympy.symbols('a D μ_0 I_0 σ_1 σ_2')

    print("The x-directed force on the conducting material in the interval s < x < 2s is:")

    # Construct and print the final equation based on the derivation that matches answer A.
    # The structure is F_x = -a*D * (mu_0/2) * (I_0^2/D^2) * (sigma_2/(sigma_1 + sigma_2))^2
    term1_str = f"-{a}*{D}"
    term2_str = f"({mu_0}/2)"
    term3_str = f"({I_0}**2/{D}**2)"
    term4_str = f"(({sigma_2})/({sigma_1} + {sigma_2}))**2"
    
    print(f"F_x = {term1_str} * {term2_str} * {term3_str} * {term4_str}")

print_force_equation()