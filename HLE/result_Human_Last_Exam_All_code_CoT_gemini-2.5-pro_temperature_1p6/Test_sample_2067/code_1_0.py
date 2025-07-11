import numpy as np

def f(t):
    """ The function f(t) as defined in the problem. """
    return np.cos(np.pi * t)**2

def h(x):
    """ y = f(sin(pi*x)) """
    return f(np.sin(np.pi * x))

def g(y):
    """ x = f(cos(2*pi*y)) """
    return f(np.cos(2 * np.pi * y))

def F(x):
    """ The composite function F(x) = g(h(x)). Solutions are fixed points of F. """
    return g(h(x))

def solve():
    """
    Solves the problem by finding the number of solutions and the number of pairs with an integer.
    """
    # Part 1: Find the total number of solutions
    # We look for roots of F(x) - x = 0 by checking for sign changes.
    x_coords = np.linspace(0, 1, 20001)
    y_coords = F(x_coords)
    
    # The number of solutions is the number of times F(x) - x crosses zero.
    # Add one for the known solution at x=1 where F(1)-1=0 without a sign change.
    diff = y_coords - x_coords
    # Find where the sign of the difference changes
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    num_solutions = len(sign_changes)

    # At x=1, F(1)=1, so F(1)-1=0. np.diff doesn't catch this if the function doesn't cross.
    # Let's check explicitly.
    if abs(F(1) - 1) < 1e-9 and not any(abs(x_coords[i] - 1) < 1e-9 for i in sign_changes):
        # We search near 1 if it is a root.
        # Check if x=1 is already counted
        is_one_counted = False
        for i in sign_changes:
            if abs(x_coords[i]-1.0)<1e-4 or abs(x_coords[i+1]-1.0)<1e-4:
                 is_one_counted = True
                 break
        if not is_one_counted:
            num_solutions +=1

    # Part 2: Find the number of pairs with at least one integer
    num_integer_pairs = 0
    checked_pairs = set()

    # Case y=0
    # From y=f(sin(pi*x)), 0 = cos(pi*sin(pi*x))^2 -> sin(pi*x)=1/2
    # x = 1/6 or x = 5/6
    # From x=f(cos(2*pi*y)), x = f(cos(0)) = f(1) = 1.
    # 1/6 != 1 and 5/6 != 1. No solution with y=0.

    # Case y=1
    # From y=f(sin(pi*x)), 1 = cos(pi*sin(pi*x))^2 -> sin(pi*x) is 0 or 1
    # sin(pi*x)=0 -> x=0,1. sin(pi*x)=1 -> x=1/2.
    # Potential x values are {0, 1/2, 1}.
    # From x=f(cos(2*pi*y)), x=f(cos(2pi))=f(1)=1.
    # So we must have x=1. The only candidate pair is (1,1).
    # Check (1,1):
    # y=f(sin(pi*1))=f(0)=1. Correct.
    # x=f(cos(2*pi*1))=f(1)=1. Correct.
    # So (1,1) is a solution.
    if (1,1) not in checked_pairs:
        num_integer_pairs += 1
        checked_pairs.add((1,1))

    # Case x=0
    # From x=f(cos(2*pi*y)), 0=... -> cos(2pi*y) = +-1/2
    # y = {1/6, 1/3, 2/3, 5/6}
    # From y=f(sin(pi*x)), y=f(sin(0))=f(0)=1.
    # y must be 1, but we found y must be in {1/6, ...}. No solution.
    
    # Case x=1 was already handled by the y=1 case, which yielded (1,1).
    
    # We found one integer pair (1,1).
    
    print(f"{num_solutions},{num_integer_pairs}")

solve()