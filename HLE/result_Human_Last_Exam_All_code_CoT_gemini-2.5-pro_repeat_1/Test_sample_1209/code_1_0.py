import sympy

def solve_relativistic_projectile():
    """
    This function uses the sympy library to derive and display the formula for the
    horizontal distance (D) traveled by a particle launched horizontally from a height (h)
    with a relativistic velocity (v0).
    """
    # Define the symbols used in the physics equations.
    # h: initial height of the cliff
    # v0: initial horizontal velocity
    # g: acceleration due to gravity
    # c: speed of light
    h, v0, g, c = sympy.symbols('h v_0 g c', positive=True, real=True)

    # Step 1: Define the initial Lorentz factor, gamma_0.
    # This factor is crucial for all relativistic calculations.
    gamma_0 = 1 / sympy.sqrt(1 - v0**2 / c**2)

    # Step 2: Define the total time of flight, T.
    # This formula is derived using conservation of energy and momentum.
    # Note: The particle's mass 'm' cancels out and does not appear in the final result.
    T = sympy.sqrt((2 * h * gamma_0) / g + h**2 / c**2)

    # Step 3: Define the horizontal distance, D.
    # This is found by integrating the horizontal velocity v_x(t) from t=0 to T.
    # The result involves the inverse hyperbolic sine function, asinh.
    D = (gamma_0 * v0 * c / g) * sympy.asinh(g * T / (c * gamma_0))

    # Display the final expression for D.
    # sympy.pprint provides a readable, formatted output for the symbolic formula.
    print("The final formula for the horizontal distance D is:")
    sympy.pprint(D)
    print("\nWhere the terms are defined as:")
    print("\nTime of flight T:")
    sympy.pprint(T)
    print("\nInitial Lorentz factor gamma_0:")
    sympy.pprint(gamma_0)

solve_relativistic_projectile()