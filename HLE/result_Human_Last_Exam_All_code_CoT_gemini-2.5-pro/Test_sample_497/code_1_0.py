import sympy

def solve_capacitance():
    """
    Calculates the value of capacitor x such that the equivalent capacitance
    of the ladder network is independent of the number of cells.
    """
    # Define the symbols for the known capacitance 'c' and the unknown 'x'.
    # Both must be positive.
    c, x = sympy.symbols('c x', positive=True)

    # The condition that the total capacitance is independent of the number of cells
    # means that the terminating capacitance 'x' must be equal to the
    # characteristic capacitance of the infinite ladder network.

    # The characteristic capacitance (which we set to x) is defined by the
    # fixed-point equation C_in = C_L = x, which leads to a quadratic equation.
    # C_in = c * (c + C_L) / (3*c + 2*C_L)
    # x = c * (c + x) / (3*c + 2*x)
    # x * (3*c + 2*x) = c^2 + c*x
    # 3*c*x + 2*x^2 = c^2 + c*x
    # 2*x^2 + 2*c*x - c^2 = 0

    print("To find the value of x, we solve the following quadratic equation for x:")
    # Printing each coefficient as requested.
    print(f"{2} * x**2 + {2}*c*x - {1}*c**2 = 0\n")

    # Define the equation using sympy
    equation = sympy.Eq(2*x**2 + 2*c*x - c**2, 0)

    # Solve the quadratic equation for x
    solutions = sympy.solve(equation, x)

    # Since capacitance must be positive, we select the positive solution.
    # The solutions are typically [-c/2 - sqrt(3)*c/2, -c/2 + sqrt(3)*c/2].
    # The second one is positive because c and sqrt(3) are positive.
    positive_solution = solutions[1]

    print("Solving for x and taking the positive solution, we get:")
    # The `sympy.pretty` function can make the output look nicer in a terminal
    # that supports unicode, but for maximum compatibility, we'll print the standard form.
    final_expression = sympy.pretty(positive_solution, use_unicode=False)
    print(f"x = {final_expression}")

solve_capacitance()