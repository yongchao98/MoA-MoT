import sympy

def solve_capacitance_problem():
    """
    This function solves for the value of capacitor x such that the equivalent
    capacitance of the ladder network is independent of the number of cells.
    """
    # Define the symbols for the characteristic capacitance (which is x) and c.
    # We define them as positive since capacitance cannot be negative.
    x, c = sympy.symbols('x c', positive=True)

    # The problem reduces to finding the fixed point of the ladder's recurrence relation,
    # which leads to the following quadratic equation for x:
    # 2*x**2 + 2*c*x - c**2 = 0
    equation = sympy.Eq(2*x**2 + 2*c*x - c**2, 0)

    # Solve the quadratic equation for x
    solutions = sympy.solve(equation, x)

    # Since we defined x and c as positive, sympy will automatically return
    # the single positive solution.
    result = solutions[0]

    print("The condition that the equivalent capacitance is independent of the number of cells")
    print("requires the terminating capacitor 'x' to be equal to the characteristic capacitance of the ladder.")
    print("\nThis leads to the quadratic equation for x: 2*x**2 + 2*c*x - c**2 = 0")
    print("\nSolving this equation for x gives the following expression:")

    # We want to print the final equation with its constituent numbers clearly shown,
    # as requested. The solution is x = c * (sqrt(3) - 1) / 2.
    sqrt_3_val = 3
    minus_1_val = 1
    denominator_val = 2

    # sympy.pretty_print(sympy.Eq(x, result))
    print(f"\nx = c * (sqrt({sqrt_3_val}) - {minus_1_val}) / {denominator_val}")


solve_capacitance_problem()