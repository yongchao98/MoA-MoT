import sympy

def solve_separatrix():
    """
    This function prints the derived equation of the separatrix for the given system of ODEs.
    The derivation steps are outlined in the text explanation.
    The separatrix is found to be the unstable manifold of the saddle point (1, -1),
    and its equation is d = -u^2.
    """
    
    # The equation of the separatrix is d = -u^2.
    # We will display this equation, explicitly showing the numbers involved as requested.
    # The equation can be written as d = -1 * u^2.
    
    coefficient = -1
    variable = 'u'
    power = 2
    
    # We use sympy to format the equation nicely.
    u = sympy.Symbol(variable)
    d = sympy.Symbol('d')
    
    equation = sympy.Eq(d, coefficient * u**power)
    
    # Print the equation with explicit numbers as requested
    print("The equation for the separatrix is:")
    print(f"d = {coefficient} * {variable}**{power}")
    
    # For a more standard mathematical output:
    print("\nIn standard mathematical notation:")
    print(sympy.pretty(equation))

solve_separatrix()