import sympy as sp

def solve_emi_shielding_force():
    """
    This function prints the symbolic formula for the force per unit area on the conductor.
    The formula is derived from first principles of magnetostatics.
    """
    # Define symbolic variables for the parameters in the problem
    mu_0, mu, K_0, a, d, y = sp.symbols('mu_0 mu K_0 a d y', real=True, positive=True)

    # The final expression for the force per unit area is derived as:
    # f/area = (mu_0 / 2) * (K_0**2 * sin(a*y)**2) / (cosh(a*d) + (mu_0 / mu) * sinh(a*d))**2 * i_x

    print("The final equation for the force per unit area on the x = d interface is:")
    print("f/area = (Coefficient) * (Numerator Term) / (Denominator Term) * (Direction Vector)\n")

    # Breaking down the final equation into its components as requested
    print("--- Components of the Final Equation ---")

    # The constant coefficient part of the equation
    coefficient = mu_0 / 2
    print(f"Coefficient: {coefficient}")

    # The numerator part of the equation, which depends on the current and position y
    numerator_term = K_0**2 * sp.sin(a*y)**2
    print(f"Numerator Term: {numerator_term}")

    # The denominator part of the equation, which depends on the material properties and geometry
    denominator_term = (sp.cosh(a*d) + (mu_0 / mu) * sp.sinh(a*d))**2
    print(f"Denominator Term: {denominator_term}")

    # The direction of the force vector
    direction_vector = "i_x (unit vector in the positive x-direction)"
    print(f"Direction Vector: {direction_vector}")
    print("----------------------------------------\n")

    # Constructing and printing the full formula for the x-component of the force
    force_x_component = coefficient * numerator_term / denominator_term
    print("The full expression for the x-component of the force per unit area is:")
    sp.pprint(force_x_component, use_unicode=True)

if __name__ == '__main__':
    solve_emi_shielding_force()