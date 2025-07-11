import numpy as np
from scipy.optimize import root

def solve_equations():
    """
    Solves the system of equations to find the size of set S and the number of pairs with an integer.
    """
    # Define the function f(t) from the problem description
    def f(t):
        return np.cos(np.pi * t)**2

    # Define the system of equations for the root finder
    # We want to find (x, y) such that F(x, y) = [0, 0]
    def equations(vars):
        x, y = vars
        eq1 = y - f(np.sin(np.pi * x))
        eq2 = x - f(np.cos(2 * np.pi * y))
        return [eq1, eq2]

    # Use a set to store unique solutions found
    solutions = set()

    # Create a grid of starting points to find all solutions in [0,1]x[0,1]
    # A 15x15 grid provides a good balance between coverage and computation time.
    grid_points = np.linspace(0, 1, 15)

    for x0 in grid_points:
        for y0 in grid_points:
            # Run the numerical solver
            sol = root(equations, [x0, y0], method='hybr', tol=1e-8)
            
            # If the solver converged successfully...
            if sol.success:
                x, y = sol.x
                # ...and the solution is within the domain [0,1]x[0,1]
                if -1e-8 <= x <= 1 + 1e-8 and -1e-8 <= y <= 1 + 1e-8:
                    # Round to 6 decimal places to group close solutions
                    solution_tuple = (round(x, 6), round(y, 6))
                    solutions.add(solution_tuple)

    # Convert the set to a sorted list for consistent output
    sorted_solutions = sorted(list(solutions))

    # The size of S is the number of unique solutions
    size_of_S = len(sorted_solutions)

    # Count how many pairs contain at least one integer coordinate (0 or 1)
    count_integer_pairs = 0
    for x, y in sorted_solutions:
        # Check if x or y is close to an integer 0 or 1
        is_x_int = np.isclose(x, 0) or np.isclose(x, 1)
        is_y_int = np.isclose(y, 0) or np.isclose(y, 1)
        if is_x_int or is_y_int:
            count_integer_pairs += 1

    # Print the results as a comma-separated list
    print(f"{size_of_S},{count_integer_pairs}")

solve_equations()