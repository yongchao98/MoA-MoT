import sympy

def solve_force_equation():
    """
    This function constructs and prints the symbolic expression for the instantaneous force
    f_x(t) based on the provided problem description and choice B.
    """
    # Define all symbolic variables
    R, N, N_0, I_0, i_0, g, mu_0, alpha_T, T, T_0, B_s, omega, t = sympy.symbols(
        'R N N_0 I_0 i_0 g mu_0 alpha_T T T_0 B_s omega t'
    )

    # Time-varying current
    i_t = i_0 * sympy.sin(omega * t)

    # Numerator of the main fraction
    # It consists of the magnetic and electrical parameters
    numerator = mu_0 * (1 - alpha_T * (T - T_0)) * N_0 * I_0 * i_t

    # Denominator of the main fraction
    # It includes the gap 'g' and the saturation effect
    saturation_term = 1 + (mu_0 * N_0 * I_0) / (g * B_s)
    denominator = g**2 * saturation_term

    # Prefactor for the force equation
    prefactor = -2 * sympy.pi * R * N

    # The full expression for the instantaneous force f_x(t)
    f_x_t = prefactor * (numerator / denominator)

    # Print the equation in a readable format
    print("The instantaneous force f_x(t) is given by the equation:")
    print("\nf_x(t) = -2*\u03C0*R*N * (numerator / denominator)\n")
    print("Where:")
    print(f"numerator = {numerator}")
    print(f"denominator = {denominator}\n")
    print("The full equation is:")
    sympy.pprint(f_x_t, use_unicode=True)

solve_force_equation()