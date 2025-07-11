import sympy

def display_fabry_perot_coefficients():
    """
    This function defines and prints the expressions for the overall reflection (rho)
    and transmission (tau) coefficients of a Fabry-Pérot cavity.
    """
    # Define the symbolic variables
    rho_m, tau_m, k0, d = sympy.symbols('ρ_m τ_m k_0 d')
    I = sympy.I

    # Phase term for a single pass
    phi_single = k0 * d
    # Phase term for a round trip
    phi_round_trip = 2 * k0 * d

    # Expression for the transmission coefficient (tau)
    numerator_tau = tau_m**2 * sympy.exp(I * phi_single)
    denominator_tau = 1 - rho_m**2 * sympy.exp(I * phi_round_trip)
    tau = numerator_tau / denominator_tau

    # Expression for the reflection coefficient (rho) from option D
    # This form is taken from the problem's chosen answer.
    # It simplifies to numerator_rho_D / denominator_rho_D
    numerator_rho_D = 1 - (rho_m - tau_m**2) * sympy.exp(I * phi_round_trip) * rho_m
    denominator_rho_D = 1 - rho_m**2 * sympy.exp(I * phi_round_trip)
    rho = numerator_rho_D / denominator_rho_D

    # Print the results in a readable format
    print("Fabry-Pérot Cavity Coefficients:")
    print("-" * 35)

    # Print the transmission coefficient equation
    print("Transmission coefficient (τ):")
    print("τ = (τ_m**2 * exp(i*k_0*d)) / (1 - ρ_m**2 * exp(i*2*k_0*d))")
    
    print("\nWhich corresponds to:")
    tau_latex = sympy.latex(tau).replace('k_{0}', 'k_0')
    # Manually format the output to match the desired style
    num_tau_str = "τ_m**2 * e**(i*k_0*d)"
    den_tau_str = "1 - ρ_m**2 * e**(i*2*k_0*d)"
    print(f"τ = ({num_tau_str}) / ({den_tau_str})")


    print("\n" + "-" * 35)

    # Print the reflection coefficient equation
    print("Reflection coefficient (ρ):")
    print("ρ = (1 - (ρ_m - τ_m**2) * exp(i*2*k_0*d) * ρ_m) / (1 - ρ_m**2 * exp(i*2*k_0*d))")

    print("\nWhich corresponds to:")
    rho_latex = sympy.latex(rho).replace('k_{0}', 'k_0')
    # Manually format the output to match the desired style
    num_rho_str = "1 - (ρ_m - τ_m**2) * e**(i*2*k_0*d) * ρ_m"
    den_rho_str = "1 - ρ_m**2 * e**(i*2*k_0*d)"
    print(f"ρ = ({num_rho_str}) / ({den_rho_str})")
    print("-" * 35)


if __name__ == '__main__':
    display_fabry_perot_coefficients()