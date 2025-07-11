import sympy

def solve_capacitor_value():
    """
    This script finds the value of capacitor x for a ladder circuit
    such that the total equivalent capacitance is independent of the number of cells.
    
    The condition for independence requires that the terminating capacitance x
    be equal to the characteristic capacitance of the infinite ladder. This leads
    to a quadratic equation for x in terms of c.
    """

    # Define symbols for the cell capacitance 'c' and the terminating capacitance 'x'.
    # Both must be positive real values.
    c, x = sympy.symbols('c x', positive=True)

    # The recurrence relation for the equivalent capacitance C_in of a cell
    # loaded by C_load is: C_in = c * (c + C_load) / (3*c + 2*C_load).
    # For independence, C_in = C_load = x.
    # This leads to the equation: x = c*(c+x) / (3c+2x)
    # which simplifies to: 2*x**2 + 2*c*x - c**2 = 0.
    
    print("The derived equation for x is:")
    equation = sympy.Eq(2*x**2 + 2*c*x - c**2, 0)
    sympy.pprint(equation, use_unicode=True)
    
    # Solve the quadratic equation for x
    solutions = sympy.solve(equation, x)
    
    print("\nSolving the equation for x gives two solutions:")
    sympy.pprint(solutions, use_unicode=True)

    # Since capacitance must be a positive quantity, we select the positive solution.
    positive_solution = None
    for sol in solutions:
        # Check if the expression is positive by substituting a positive value for c (e.g., c=1)
        if sol.subs(c, 1) > 0:
            positive_solution = sol
            break
    
    print("\nThe physically valid (positive) solution for x is:")
    sympy.pprint(positive_solution, use_unicode=True)

    # The prompt requests to output each number in the final equation.
    # The solution is of the form: x = c * (sqrt(A) - B) / D
    # For our solution x = c*(-1 + sqrt(3))/2, the values are:
    A = 3
    B = 1
    D = 2
    
    print("\nThe final equation can be written as: x = c * (sqrt(A) - B) / D")
    print("The numbers in this final equation are:")
    print(f"A (inside the square root) = {A}")
    print(f"B (subtracted from the square root) = {B}")
    print(f"D (in the denominator) = {D}")

if __name__ == '__main__':
    solve_capacitor_value()