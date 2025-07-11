import sympy

def solve_h_function():
    """
    This function derives and prints the expression for h(x).
    """
    # Define x as a symbolic variable representing b(0)
    x = sympy.symbols('x')

    # The condition for a(t) to converge to 0 is that a'(0) > 0.
    # a'(0) = -1/2*a(0)**2 - 2*b(0)**2 - b(0) + 1 > 0
    # Rearranging for a(0)**2 gives:
    # a(0)**2 < 2 - 2*b(0) - 4*b(0)**2
    #
    # The problem asks for a condition of the form -sqrt(h(b(0))) < a(0) < 0.
    # For a(0) < 0, this is equivalent to a(0)**2 < h(b(0)).
    # So, we identify h(x) where x = b(0).
    
    constant_term = 2
    linear_coeff = -2
    quadratic_coeff = -4
    
    h_x = constant_term + linear_coeff*x + quadratic_coeff*x**2
    
    print("The function h(x) is derived from the condition that a'(0) must be positive.")
    print("This defines a region in the phase space from which trajectories converge to a(t) -> 0.")
    print("\nThe derived function h(x) is:")
    
    # We output each number and the final equation as requested.
    equation_string = f"h(x) = {constant_term} + ({linear_coeff})*x + ({quadratic_coeff})*x**2"
    
    print(equation_string)

solve_h_function()