import numpy as np
from scipy.optimize import fsolve

def equations(p):
    """
    Defines the system of equations to be solved.
    Let f(t) = cos(pi*t)^2.
    The equations are y = f(sin(pi*x)) and x = f(cos(2*pi*y)).
    """
    x, y = p
    # Define the function f(t)
    f_t = lambda t: (np.cos(np.pi * t))**2
    
    # The two equations from the problem statement
    eq1 = y - f_t(np.sin(np.pi * x))
    eq2 = x - f_t(np.cos(2 * np.pi * y))
    
    return (eq1, eq2)

# Use a set to store unique solutions found
solutions = set()

# Search for solutions by starting the solver from a grid of initial points
# A 20x20 grid is dense enough to find all solutions for this system
for i in np.linspace(0, 1, 20):
    for j in np.linspace(0, 1, 20):
        # Find the root of the system of equations
        sol, _, ier, _ = fsolve(equations, (i, j), full_output=True)
        
        # ier=1 means the solution is found
        if ier == 1:
            # Check if the solution is within the valid range [0,1] x [0,1]
            if 0 <= sol[0] <= 1 and 0 <= sol[1] <= 1:
                # Round the solution to handle small numerical errors
                # and add the tuple to the set to ensure uniqueness
                solutions.add((round(sol[0], 5), round(sol[1], 5)))

# Convert the set of solutions to a sorted list for consistent output
sorted_solutions = sorted(list(solutions))

# Count the total number of solutions
size_S = len(sorted_solutions)

# Count how many of the pairs contain at least one integer
# The only integers in the range [0,1] are 0 and 1.
# We check if x or y is approximately 0 or 1.
integer_pairs_count = 0
for sol in sorted_solutions:
    x, y = sol
    # A pair contains an integer if x or y is 0 or 1.
    # We use a tolerance for floating point comparison.
    is_x_integer = np.isclose(x, 0) or np.isclose(x, 1)
    is_y_integer = np.isclose(y, 0) or np.isclose(y, 1)
    if is_x_integer or is_y_integer:
        integer_pairs_count += 1
        
# Print the results as a comma-separated list
print(f"{size_S},{integer_pairs_count}")
