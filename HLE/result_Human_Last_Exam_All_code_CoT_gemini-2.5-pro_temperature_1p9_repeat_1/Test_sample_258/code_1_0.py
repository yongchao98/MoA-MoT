import math

def solve_grid_crossing():
    """
    Calculates the minimal and maximal number of grid cells a circle
    of radius 500 can cross, based on its center position.
    """
    R = 500
    
    # --- Minimal Number of Crossings ---
    # This occurs when the center (xc, yc) is placed infinitesimally close to
    # an integer grid point, e.g., (epsilon, epsilon).
    # For a given dimension, let's take xc = epsilon.
    # The number of vertical integer lines crossed is the number of integers k
    # in the open interval (xc - R, xc + R) -> (epsilon - 500, epsilon + 500).
    # The integers are -499, -498, ..., 499.
    min_n_per_dimension = 499 - (-499) + 1
    
    min_vertical_crossings = 2 * min_n_per_dimension
    min_horizontal_crossings = 2 * min_n_per_dimension
    min_total_cells = min_vertical_crossings + min_horizontal_crossings

    print("To find the minimum number of cells crossed:")
    print(f"We choose a center close to an integer coordinate (e.g., (0.001, 0.001)).")
    print(f"Number of vertical lines crossed (n_v) = {min_n_per_dimension}")
    print(f"Each line is crossed twice, so vertical crossings = 2 * {min_n_per_dimension} = {min_vertical_crossings}")
    print(f"Number of horizontal lines crossed (n_h) = {min_n_per_dimension}")
    print(f"Each line is crossed twice, so horizontal crossings = 2 * {min_n_per_dimension} = {min_horizontal_crossings}")
    print(f"Minimal number of cells = {min_vertical_crossings} + {min_horizontal_crossings} = {min_total_cells}\n")

    # --- Maximal Number of Crossings ---
    # This occurs when the center (xc, yc) is placed at a half-integer
    # coordinate, e.g., (0.5, 0.5).
    # For a given dimension, let's take xc = 0.5.
    # The number of vertical integer lines crossed is the number of integers k
    # in the open interval (xc - R, xc + R) -> (0.5 - 500, 0.5 + 500) -> (-499.5, 500.5).
    # The integers are -499, -498, ..., 500.
    max_n_per_dimension = 500 - (-499) + 1

    max_vertical_crossings = 2 * max_n_per_dimension
    max_horizontal_crossings = 2 * max_n_per_dimension
    max_total_cells = max_vertical_crossings + max_horizontal_crossings

    print("To find the maximum number of cells crossed:")
    print(f"We choose a center at a half-integer coordinate (e.g., (0.5, 0.5)).")
    print(f"Number of vertical lines crossed (n_v) = {max_n_per_dimension}")
    print(f"Each line is crossed twice, so vertical crossings = 2 * {max_n_per_dimension} = {max_vertical_crossings}")
    print(f"Number of horizontal lines crossed (n_h) = {max_n_per_dimension}")
    print(f"Each line is crossed twice, so horizontal crossings = 2 * {max_n_per_dimension} = {max_horizontal_crossings}")
    print(f"Maximal number of cells = {max_vertical_crossings} + {max_horizontal_crossings} = {max_total_cells}\n")
    
    final_answer = (min_total_cells, max_total_cells)
    print(f"The final answer in the form (min, max) is: {final_answer}")

solve_grid_crossing()