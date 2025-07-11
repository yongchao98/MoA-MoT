import sympy

def solve_capacitance():
    """
    This function solves for the capacitance x in terms of c for the given ladder network.
    """
    # Define symbolic variables for c and x. We assume c is positive.
    c = sympy.Symbol('c', positive=True)
    x = sympy.Symbol('x')

    # The condition for the equivalent capacitance to be independent of the number of cells
    # is that x is the characteristic capacitance of the ladder. This leads to the
    # quadratic equation for x: 8*x**2 + 3*c*x - 2*c**2 = 0.
    
    # We can solve this by defining a ratio r = x/c, which gives 8*r**2 + 3*r - 2 = 0.
    # The coefficients for the quadratic formula a*r**2 + b*r + c = 0 are:
    a = 8
    b = 3
    c_coeff = -2
    
    print("The problem reduces to finding the characteristic capacitance 'x' of the ladder.")
    print("This leads to the following quadratic equation for x in terms of c:")
    # Using sympy.pretty to print the equation
    equation = sympy.Eq(8 * x**2 + 3 * c * x - 2 * c**2, 0)
    print("Equation:", sympy.pretty(equation))
    
    print("\nTo solve this, we can find the ratio r = x/c.")
    print(f"The equation for r is: {a}*r^2 + {b}*r + {c_coeff} = 0")
    
    print("\nUsing the quadratic formula r = (-b +/- sqrt(b^2 - 4*a*c)) / (2*a):")
    
    discriminant = b**2 - 4 * a * c_coeff
    
    # Printing the substitution of numbers into the formula
    print(f"r = ( -({b}) +/- sqrt( ({b})^2 - 4*({a})*({c_coeff}) ) ) / ( 2*({a}) )")
    print(f"r = ( -{b} +/- sqrt( {b**2} - ({4*a*c_coeff}) ) ) / {2*a}")
    print(f"r = ( -{b} +/- sqrt( {discriminant} ) ) / {2*a}")

    # Since capacitance must be positive, x must be positive, and thus r must be positive.
    # We take the positive root.
    r_sol_pos = (-b + sympy.sqrt(discriminant)) / (2*a)
    
    print("\nSince capacitance must be positive, we take the positive root for r:")
    print(f"r = ( -{b} + sqrt({discriminant}) ) / {2*a}")
    
    # Final symbolic solution for x
    x_solution = c * r_sol_pos
    
    print("\nThe final value for x in terms of c is:")
    # Using sympy's pretty print for the final expression
    print(sympy.pretty(sympy.Eq(x, x_solution)))
    
solve_capacitance()