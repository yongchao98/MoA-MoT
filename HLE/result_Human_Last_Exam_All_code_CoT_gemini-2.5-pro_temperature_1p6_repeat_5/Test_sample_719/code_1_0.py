import sympy

def solve_problem():
    """
    This script defines the mathematical objects as symbolic variables
    and prints the derived expression for theta'(t).
    """
    # Define t as a symbol
    t = sympy.Symbol('t')
    
    # Define f and theta as functions of t
    f = sympy.Function('f')(t)
    theta = sympy.Function('theta')(t)
    
    # The first term of the final equation for theta'(t)
    term1 = f * sympy.cos(theta)**2
    
    # The second term of the final equation for theta'(t)
    term2 = sympy.diff(f, t)/f * sympy.cos(theta) * sympy.sin(theta)
    
    # The full expression for theta'(t)
    theta_prime = term1 + term2
    
    print("The derived expression for theta'(t) is:")
    sympy.pprint(theta_prime)
    
    print("\nLet's print the terms of the final equation separately:")
    
    print("\nTerm 1:")
    sympy.pprint(term1)

    print("\nTerm 2:")
    sympy.pprint(term2)
    
solve_problem()