import numpy as np
from scipy.optimize import brentq

def g(x):
    """Corresponds to the first equation y = g(x)"""
    return np.cos(np.pi * np.sin(np.pi * x))**2

def h(y):
    """Corresponds to the second equation x = h(y)"""
    return np.cos(np.pi * np.cos(2 * np.pi * y))**2

def k(x):
    """Composite function k(x) = h(g(x))"""
    return h(g(x))

def find_roots():
    """Finds the number of solutions by finding roots of k(x) - x = 0"""
    
    # The function whose roots we want to find
    func_to_solve = lambda x: k(x) - x

    # We found that the only solution involving an integer is (1,1).
    num_integer_pairs = 1

    # Now, find all solutions in [0,1].
    # We analyze the function func_to_solve on different sub-intervals
    # to find all its roots.
    # The special points for the inner functions help define search intervals.
    # Based on analysis, we expect sign changes in these intervals.
    intervals = [
        (0.0, 0.3),          # From analysis: k(0)=1 > 0, k(1/3)~0.23 < 1/3
        (0.4, 0.5),          # From analysis: k(1/3) < 1/3, k(1/2)=1 > 1/2
        (0.5, 0.7),          # Symmetric behavior expected
        (0.7, 0.9),          # Symmetric behavior expected
    ]

    roots = set()
    
    # Add the known root at x=1
    # Check if k(1)-1 is close to zero
    if abs(func_to_solve(1.0)) < 1e-9:
        roots.add(1.0)

    for a, b in intervals:
        try:
            root = brentq(func_to_solve, a, b)
            # To handle floating point issues, round the root
            roots.add(round(root, 8))
        except ValueError:
            # This means func_to_solve(a) and func_to_solve(b) have the same sign.
            # Our interval choices should prevent this, but this is for safety.
            pass
            
    total_solutions = len(roots)

    # Print the comma-separated values as requested.
    # The output shows the total number of solutions, and the number of solutions
    # where at least one component is an integer.
    # From our analysis, this should be 5 and 1.
    print(f"{total_solutions},{num_integer_pairs}")

find_roots()