import sympy

def solve_mass_equation():
    """
    This function symbolically calculates the squared mass of the sixth degree of freedom.
    """
    # Define the symbols we will use in our equations.
    # M_sq will represent the squared mass we are looking for.
    # m is the mass parameter from the Lagrangian.
    M_sq, m = sympy.symbols('M_sq m')
    
    # h represents the scalar field (trace of the metric perturbation).
    # Box is a placeholder for the d'Alembert wave operator.
    h = sympy.Function('h')
    Box = sympy.Function('Box')

    # From the theoretical derivation, the equation of motion for the scalar mode h is:
    # 2 * Box(h) + m**2 * h = 0
    # We represent this equation symbolically.
    scalar_eom = sympy.Eq(2 * Box(h) + m**2 * h, 0)

    # The standard Klein-Gordon equation for a particle with squared mass M_sq is:
    # Box(h) + M_sq * h = 0
    # Our goal is to rearrange our derived EOM into this standard form to identify M_sq.

    # We can solve for Box(h) from our EOM
    # 2*Box(h) = -m**2 * h
    # Box(h) = -(m**2 / 2) * h
    # Comparing this to Box(h) = -M_sq * h, we find M_sq.
    
    # We can use sympy to solve for the term that corresponds to M_sq
    # From Box(h) + (m**2 / 2) * h = 0
    result_squared_mass = (m**2 / 2)

    # Now, we print the result clearly, showing each number in the final equation.
    # The result is of the form C * m**p
    coeff = result_squared_mass.as_coeff_mul()[0] # This gives the fraction 1/2
    power = sympy.degree(result_squared_mass, m) # This gives the power of m, which is 2
    
    numerator = coeff.p
    denominator = coeff.q
    
    print("The theoretical derivation leads to the equation of motion for the scalar field h (the 6th degree of freedom):")
    print(f"    2 * Box(h) + m^2 * h = 0")
    print("\nTo find the squared mass, we compare this to the standard Klein-Gordon equation, Box(h) + M_sq * h = 0.")
    print("Rearranging our equation gives:")
    print(f"    Box(h) + ({numerator}/{denominator}) * m^{power} * h = 0")
    print("\nTherefore, the squared mass of the sixth degree of freedom is:")
    print(f"M_sq = {numerator}/{denominator} * m^{power}")

# Run the solver
solve_mass_equation()