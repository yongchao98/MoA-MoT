import sympy

def solve_ladder_capacitance():
    """
    Solves for the terminating capacitance x that makes the equivalent
    capacitance of a ladder circuit independent of the number of cells.
    """
    # Step 1: Explain the problem setup and the condition for the solution.
    print("This script solves for the terminating capacitance 'x' for the given ladder circuit.")
    print("The goal is to find a value for x (in terms of c) such that the total")
    print("equivalent capacitance between A and B is independent of the number of cells N.")
    print("\nThis condition is met if, when we add one cell to a ladder with equivalent")
    print("capacitance x, the new equivalent capacitance is also x.")
    print("-" * 60)

    # Step 2: Define symbols and the governing equation.
    c, x = sympy.symbols('c x', positive=True)

    print("After analyzing the circuit of a single cell terminated by a capacitor 'x',")
    print("we find the recurrence relation for the new equivalent capacitance (C_new):")
    print("C_new = (c * (c + x)) / (3*c + 2*x)")

    print("\nTo make the total capacitance independent of the number of cells, we set C_new = x:")
    # Print the equation using the symbols
    print(f"       x = (c * (c + x)) / (3*c + 2*x)")

    # Step 3: Rearrange into a quadratic equation and show it.
    print("\nRearranging this gives the following quadratic equation for x:")
    # The numbers in this equation are 2, 2, and -1.
    print(f"       2*x**2 + 2*c*x - 1*c**2 = 0")

    # Step 4: Solve the equation and display the result.
    # Set up the equation in SymPy
    equation = sympy.Eq(2*x**2 + 2*c*x - c**2, 0)
    
    # Solve for x
    solutions = sympy.solve(equation, x)
    
    # SymPy automatically uses the 'positive=True' assumption to return only the valid solution.
    final_solution = solutions[0]

    print("\nSolving this quadratic equation for x (and taking the positive root, as capacitance cannot be negative) yields the final answer.")
    print("\nThe numbers in the final equation for x are 1, 3, and 2.")
    print("\nThe required value for the capacitor x is:")
    
    # We want to print the answer in a clear, human-readable format.
    # The solution from sympy is c*(sqrt(3)/2 - 1/2)
    print("       x = c * (sqrt(3) - 1) / 2")
    
    # We can also print the nicely formatted SymPy object for verification
    # print("\nFormatted SymPy output:")
    # sympy.pprint(final_solution, use_unicode=True)

# Run the solver function
solve_ladder_capacitance()
