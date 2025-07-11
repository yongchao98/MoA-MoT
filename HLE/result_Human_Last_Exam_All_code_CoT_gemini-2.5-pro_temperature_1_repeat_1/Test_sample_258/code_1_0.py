import math

def solve_grid_crossing():
    """
    Calculates the minimal and maximal numbers of grid cells a circle
    of radius 500 can cross based on established geometric formulas.
    """
    R = 500

    # For a large integer radius R, the maximal number of cells crossed
    # by a circle is given by the formula 8R. This typically occurs when
    # the circle's center is infinitesimally close to a grid point.
    max_cells = 8 * R
    
    # The minimal number of cells crossed is given by the formula 8R - 8.
    # This typically occurs when the circle's center is at (k+0.5, l+0.5),
    # maximizing symmetry and path overlap.
    min_cells = 8 * R - 8
    
    print(f"Radius (R) = {R}")
    print(f"Maximal number of cells crossed = 8 * R = 8 * {R} = {max_cells}")
    print(f"Minimal number of cells crossed = 8 * R - 8 = 8 * {R} - 8 = {min_cells}")
    print(f"The answer is ({min_cells}, {max_cells}).")

solve_grid_crossing()
