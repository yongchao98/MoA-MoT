import math

def solve_grid_circle_problem():
    """
    Calculates the minimal and maximal numbers of grid cells a circle of radius 500 can cross.
    """
    R = 500

    # Calculate the maximal number of cells
    # This occurs when the center is infinitesimally close to a grid node.
    # The formula is N_max = 8R - 8 for integer R > 1.
    max_val = 8 * R - 8
    print(f"For the maximal case:")
    print(f"The number of cells is given by the formula 8 * R - 8.")
    print(f"N_max = 8 * {R} - 8 = {max_val}")
    print("-" * 20)

    # Calculate the minimal number of cells
    # This occurs when the center is placed to maximize corner-cutting.
    # The formula is N_min = 8R - 4 * floor((sqrt(2) - 1)R) - 4.
    sqrt2_minus_1_R = (math.sqrt(2) - 1) * R
    floor_val = math.floor(sqrt2_minus_1_R)
    min_val = 8 * R - 4 * floor_val - 4
    
    print(f"For the minimal case:")
    print(f"The number of cells is given by the formula 8 * R - 4 * floor((sqrt(2) - 1) * R) - 4.")
    print(f"N_min = 8 * {R} - 4 * floor(({math.sqrt(2):.6f} - 1) * {R}) - 4")
    print(f"N_min = 8 * {R} - 4 * floor({sqrt2_minus_1_R:.6f}) - 4")
    print(f"N_min = 8 * {R} - 4 * {floor_val} - 4")
    print(f"N_min = {8 * R} - {4 * floor_val} - 4 = {min_val}")
    print("-" * 20)

    print(f"The minimal and maximal numbers of grid's cells are ({min_val}, {max_val}).")

solve_grid_circle_problem()