import sympy

def display_force_equation():
    """
    This function defines the variables symbolically and prints the equation for
    the instantaneous force fx(t) based on the derived formula.
    """
    # Define all variables as symbolic objects
    t, R, N, N_0, I_0, i_0, omega = sympy.symbols('t R N N_0 I_0 i_0 omega', real=True, positive=True)
    g, mu_0, alpha_T, T, T_0, B_s = sympy.symbols('g mu_0 alpha_T T T_0 B_s', real=True, positive=True)

    # Time-varying current i(t)
    i_t = i_0 * sympy.sin(omega * t)

    # Temperature correction term for permeability
    temp_term = 1 - alpha_T * (T - T_0)

    # Saturation term in the denominator
    saturation_term = 1 + (mu_0 * N_0 * I_0) / (g * B_s)

    # Numerator of the force equation
    numerator = -2 * sympy.pi * R * N * mu_0 * temp_term * N_0 * I_0 * i_t

    # Denominator of the force equation
    denominator = g**2 * saturation_term

    # Full expression for the instantaneous force fx(t)
    fx_t = numerator / denominator

    # Print the equation in a readable format
    print("The derived formula for the instantaneous force fx(t) is:")
    sympy.pprint(fx_t, use_unicode=True)

if __name__ == '__main__':
    display_force_equation()
