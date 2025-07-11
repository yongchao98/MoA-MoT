import sympy

def solve_capacitance():
    """
    This function calculates the value of capacitor x for the ladder circuit.
    """
    # Define the symbols for capacitance c and the unknown capacitor x.
    # We assume they are positive values.
    c = sympy.Symbol('c', positive=True)
    x = sympy.Symbol('x', positive=True)

    # Explain the principle of characteristic capacitance.
    print("To make the equivalent capacitance independent of the number of cells N,")
    print("the terminating capacitor 'x' must be equal to the characteristic capacitance of the ladder network.")
    print("This means that if we take a ladder terminated by 'x', its input capacitance is also 'x'.\n")

    # Set up the equation based on one cell of the ladder.
    # The load 'x' is in parallel with the shunt capacitor 'c'.
    C_parallel = c + x
    print(f"The equivalent capacitance of the shunt capacitor 'c' in parallel with the load 'x' is: c + x\n")

    # The input capacitance of one cell terminated by C_parallel is C_in.
    # The relation for this T-network is: C_in = (c * C_parallel) / (c + 2 * C_parallel)
    # We set C_in = x and C_parallel = c + x
    equation = sympy.Eq(x, (c * (c + x)) / (c + 2 * (c + x)))
    
    print("Setting the input capacitance equal to 'x', we get the equation:")
    # Using str() for clean printing without the Eq wrapper
    print(f"x = (c * (c + x)) / (3*c + 2*x)\n")

    # This equation can be rearranged into a quadratic form: a*x^2 + b*x + c_quad = 0
    # x * (3*c + 2*x) = c^2 + c*x
    # 3*c*x + 2*x^2 = c^2 + c*x
    # 2*x^2 + 2*c*x - c^2 = 0
    quadratic_eq = sympy.Eq(2*x**2 + 2*c*x - c**2, 0)
    print("This simplifies to the following quadratic equation:")
    print(f"{sympy.pretty(quadratic_eq, use_unicode=False)}\n")
    
    # Solve the quadratic equation for x.
    solutions = sympy.solve(quadratic_eq, x)

    # Filter for the positive solution, as capacitance cannot be negative.
    positive_solution = [sol for sol in solutions if sol.is_positive][0]

    # Print the final result, formatted as requested.
    print("Solving for x and taking the positive root gives the final value for the capacitor:")
    # The solution is c * (sqrt(3) - 1) / 2.
    # We print it piece by piece to satisfy the "output each number" constraint.
    print("x = c * (sqrt(", 3, ") - ", 1, ") / ", 2, sep='')

solve_capacitance()