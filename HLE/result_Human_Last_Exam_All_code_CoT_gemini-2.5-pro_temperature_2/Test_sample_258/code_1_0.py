import math

def solve_circle_grid_problem(radius):
    """
    Calculates the minimal and maximal number of grid cells crossed by a circle.

    Args:
        radius (int): The radius of the circle.

    Returns:
        tuple: A tuple containing the minimal and maximal number of cells.
    """
    R = int(radius)

    # For a circle of integer radius R, placed on a grid such that it's not tangent
    # to any grid lines, the number of line crossings is constant.
    # Number of vertical lines crossed = 2 * R
    # Number of horizontal lines crossed = 2 * R
    # Each line is crossed twice.
    # Total line crossings = 2 * (2*R) + 2 * (2*R) = 8 * R
    
    # The maximum number of cells is crossed when the circle's center is in a generic position
    # (e.g., (epsilon, epsilon) relative to a grid point).
    # In this case, the number of cells crossed equals the number of line crossings.
    maximal_cells = 8 * R
    
    # The minimum number of cells is crossed when the circle's center is at a point of high
    # symmetry, e.g., (k+0.5, l+0.5) for integers k, l. For an integer radius, this
    # leads to a reduction of 8 cells from the maximum. This is because at each of the
    # four extremal points of the circle, the symmetric path saves two cells.
    minimal_cells = 8 * R - 8
    
    print(f"The radius of the circle is R = {R}.")
    print(f"The maximal number of cells crossed is calculated as 8 * R.")
    print(f"8 * {R} = {maximal_cells}")
    print(f"The minimal number of cells crossed is calculated as 8 * R - 8.")
    print(f"8 * {R} - 8 = {minimal_cells}")
    
    print("\nThe final answer in the format (min, max) is:")
    print(f"({minimal_cells}, {maximal_cells})")
    
    return minimal_cells, maximal_cells

# The radius given in the problem
R = 500
solve_circle_grid_problem(R)