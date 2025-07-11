import math

def solve():
    """
    Calculates the minimal and maximal number of grid cells crossed by a circle.
    
    The number of cells a simple closed convex curve (like a circle) crosses on a grid,
    given it does not pass through any grid vertices, is N = V + H.
    V = total crossings with vertical lines.
    H = total crossings with horizontal lines.
    
    Since a circle crosses each line it intersects twice (no tangency),
    V = 2 * n_v and H = 2 * n_h, where n_v and n_h are the number of
    vertical and horizontal grid lines intersected.
    N = 2 * n_v + 2 * n_h.
    
    n_v is the number of integers in the interval [x_c - R, x_c + R].
    Since R is an integer (500) and x_c is not an integer (due to no tangency),
    the number of integers in this interval of length 2R is exactly 2R.
    So, n_v = 2R and n_h = 2R.
    
    The total number of cells is N = 2 * (2R) + 2 * (2R) = 8R.
    As this value is constant, the minimum and maximum are the same.
    """
    R = 500
    
    # n_v is the number of vertical grid lines crossed.
    # As derived in the explanation, n_v = 2 * R because R is an integer.
    n_v = 2 * R
    
    # n_h is the number of horizontal grid lines crossed.
    # Similarly, n_h = 2 * R.
    n_h = 2 * R
    
    # V is the total number of vertical crossings.
    # Each of the n_v lines is crossed twice.
    V = 2 * n_v
    
    # H is the total number of horizontal crossings.
    # Each of the n_h lines is crossed twice.
    H = 2 * n_h
    
    # The total number of cells crossed is N = V + H.
    N = V + H
    
    min_cells = N
    max_cells = N
    
    print(f"Radius (R) = {R}")
    print(f"Number of vertical lines intersected (n_v) = 2 * R = {n_v}")
    print(f"Number of horizontal lines intersected (n_h) = 2 * R = {n_h}")
    print(f"Total vertical crossings (V) = 2 * n_v = {V}")
    print(f"Total horizontal crossings (H) = 2 * n_h = {H}")
    print(f"Total cells crossed (N) = V + H = {V} + {H} = {N}")
    print(f"Since the number of cells is constant, min = max.")
    
    # Output the final answer in the required format.
    print(f"The minimal and maximal numbers are ({min_cells}, {max_cells}).")

solve()
<<< (4000, 4000) >>>