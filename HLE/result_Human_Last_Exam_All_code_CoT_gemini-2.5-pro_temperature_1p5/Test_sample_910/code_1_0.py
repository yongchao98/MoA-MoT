import sympy as sp

def solve_liquid_rise():
    """
    This function defines and prints the symbolic expressions for the liquid rise height (xi)
    and the required voltage (V0) for the height to be half the gap separation, based on Answer C.
    """
    # Define the symbolic variables
    s = sp.Symbol('s', positive=True)          # Separation between plates
    xi = sp.Symbol('xi', positive=True)         # Height of liquid rise
    epsilon_0 = sp.Symbol('varepsilon_0', positive=True) # Permittivity of free space
    V_0 = sp.Symbol('V_0', positive=True)         # Applied voltage
    rho = sp.Symbol('rho', positive=True)         # Mass density of the liquid
    g = sp.Symbol('g', positive=True)           # Gravitational acceleration
    gamma = sp.Symbol('gamma', positive=True)       # Surface tension

    # Expression for the height xi as per Answer C
    # xi = s * (epsilon_0 * V_0**2 / (2 * rho * g * s**3) - gamma / (rho * g * s))
    term1_xi = epsilon_0 * V_0**2 / (2 * rho * g * s**3)
    term2_xi = gamma / (rho * g * s)
    xi_expression = s * (term1_xi - term2_xi)

    # Expression for the voltage V_0 when xi = s/2 as per Answer C
    # V_0 = sqrt(4 * rho * g * s**3 / epsilon_0) * (1 + 2 * gamma * s / (rho * g))**(1/2)
    term1_V0 = 4 * rho * g * s**3 / epsilon_0
    term2_V0 = 1 + 2 * gamma * s / (rho * g)
    V0_expression = sp.sqrt(term1_V0 * term2_V0)


    print("Chosen Answer: C")
    print("-" * 30)
    print("An expression for the height xi of the liquid rise is:")
    sp.pprint(sp.Eq(xi, xi_expression), use_unicode=True)
    print("\nThe voltage V_0 when the liquid rise is xi = s/2 is:")
    sp.pprint(sp.Eq(V_0, V0_expression), use_unicode=True)
    print("\nStability Discussion:")
    print("The interface becomes unstable if the electrostatic forces overcome the restoring gravitational and surface tension forces. This can lead to an uncontrolled rise of the liquid, often termed pull-in instability.")

solve_liquid_rise()