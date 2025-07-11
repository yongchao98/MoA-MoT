def solve_circle_grid_problem():
    """
    Calculates the minimal and maximal number of grid cells a circle
    of radius 500 can cross based on known geometric results.
    """
    R = 500

    # The maximal number of cells crossed by a circle of integer radius R
    # is 8*R. This occurs when the center is very close to an integer coordinate.
    max_cells = 8 * R

    # The minimal number of cells is 8*R - 8. This occurs when the center
    # is at a half-integer coordinate (e.g., (0.5, 0.5)).
    min_cells = 8 * R - 8

    print(f"R = {R}")
    print(f"Maximal number of cells = 8 * {R} = {max_cells}")
    print(f"Minimal number of cells = 8 * {R} - 8 = {min_cells}")
    print(f"The minimal and maximal numbers of grid's cells are ({min_cells}, {max_cells}).")

solve_circle_grid_problem()