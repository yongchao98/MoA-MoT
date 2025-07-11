import sympy

def find_capacitor_value():
    """
    This script calculates the value of capacitor x for a ladder circuit
    such that the total equivalent capacitance is independent of the number of cells.
    """

    # Define symbolic variables. 'c' is the capacitance in the cell,
    # and 'x' is the unknown capacitance we are solving for.
    # Both must be positive.
    c = sympy.Symbol('c', positive=True)
    x = sympy.Symbol('x', positive=True)

    # The condition that the total capacitance is independent of the number
    # of cells requires the termination capacitor 'x' to be equal to the
    # characteristic capacitance of the infinite ladder network.

    # The recurrence relation for the input capacitance (C_in) of a single
    # H-network cell loaded with capacitance x is:
    # C_in = c * (c + x) / (3*c + 2*x)

    # The characteristic capacitance is the fixed point where C_in = x.
    # x = c * (c + x) / (3*c + 2*x)
    # x * (3*c + 2*x) = c**2 + c*x
    # 3*c*x + 2*x**2 = c**2 + c*x
    # This simplifies to the following quadratic equation for x:
    # 2*x**2 + 2*c*x - c**2 = 0
    
    equation = sympy.Eq(2*x**2 + 2*c*x - c**2, 0)

    # Solve the quadratic equation for x.
    solutions = sympy.solve(equation, x)

    # The physically valid solution must be positive.
    positive_solution = [sol for sol in solutions if sol.is_positive][0]

    # The final answer is the symbolic expression for x.
    # To output each number in the equation as requested, we can
    # deconstruct the symbolic expression.
    coeff = positive_solution / c
    numer, denom = sympy.fraction(coeff)
    # The numerator is of the form: sqrt(A) - B
    term1 = numer.as_ordered_terms()[0] # This will be sqrt(3)
    term2 = numer.as_ordered_terms()[1] # This will be -1

    sqrt_val = term1.args[0]
    const_val = abs(term2)

    print("The value of capacitor x must be equal to the characteristic capacitance of the ladder network.")
    print("The resulting equation for x is:")
    print(f"x = c * (sqrt({sqrt_val}) - {const_val}) / {denom}")

find_capacitor_value()