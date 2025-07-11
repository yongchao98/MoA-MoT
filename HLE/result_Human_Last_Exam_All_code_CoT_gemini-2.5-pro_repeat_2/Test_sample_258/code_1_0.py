import math

def solve_grid_crossings():
    """
    Calculates the minimal and maximal numbers of grid cells a circle of radius 500 can cross.
    """
    R = 500
    
    # The maximum number of cells crossed by a circle of radius R is 8R.
    # This occurs when the circle's center is placed at (k+0.5, l+0.5) for integers k, l.
    max_cells = 8 * R
    
    # The minimum number of cells crossed is given by 8R - 4k, where k is
    # the number of pairs of positive integers (m, n) with n > m > 0
    # such that m^2 + n^2 = R^2. This corresponds to placing the circle's center
    # infinitesimally close to a grid vertex.
    
    k = 0
    # We need to find integer solutions for m^2 + n^2 = R^2.
    # We can iterate m from 1 up to floor(R/sqrt(2)) to find unique pairs (m, n) with m < n.
    limit = int(R / math.sqrt(2))
    for m in range(1, limit + 1):
        n_squared = R*R - m*m
        n = math.isqrt(n_squared)
        
        # Check if n_squared is a perfect square and m^2 + n^2 = R^2
        if n*n == n_squared:
            # We found a Pythagorean triple (m, n, R).
            # The condition n > m is implicitly satisfied because m <= R/sqrt(2).
            k += 1
            
    min_cells = 8 * R - 4 * k
    
    print(f"The radius of the circle is R = {R}.")
    print(f"The maximum number of cells is 8 * R = 8 * {R} = {max_cells}.")
    print(f"The number of Pythagorean triples (m, n, R) with R > n > m > 0 is k = {k}.")
    print(f"The minimum number of cells is 8 * R - 4 * k = 8 * {R} - 4 * {k} = {min_cells}.")
    print(f"The minimal and maximal numbers of cells are ({min_cells}, {max_cells}).")

solve_grid_crossings()
<<< (3988, 4000) >>>