import math

def solve_grid_crossing():
    """
    Calculates the minimal and maximal number of grid cells a circle can cross.
    """
    R = 500

    # To find the maximum number of cells, we need to minimize the subtracted term.
    # This happens when the fractional parts of the center (fx, fy) are not 0.5.
    # S_fx = floor(2*fx) + floor(2*(1-fx)) = 1 for fx != 0.5
    # S_fy = floor(2*fy) + floor(2*(1-fy)) = 1 for fy != 0.5
    # We can choose a center like (0.2, 0.3)
    s_sum_min = 1 + 1
    max_cells = 8 * R - 2 * s_sum_min
    
    print(f"The radius of the circle is R = {R}.")
    print("\nTo find the maximal number of cells, we choose the center's fractional parts away from 0.5.")
    print(f"Maximal cells = 8 * R - 2 * (1 + 1) = 8 * {R} - 4")
    print(f"Maximal cells = {8*R} - 4 = {max_cells}")

    # To find the minimum number of cells, we need to maximize the subtracted term.
    # This happens when the fractional parts of the center (fx, fy) are 0.5.
    # S_fx = floor(2*0.5) + floor(2*(1-0.5)) = 1 + 1 = 2
    # S_fy = floor(2*0.5) + floor(2*(1-0.5)) = 1 + 1 = 2
    # This corresponds to a center at (I+0.5, J+0.5).
    s_sum_max = 2 + 2
    min_cells = 8 * R - 2 * s_sum_max

    print("\nTo find the minimal number of cells, we choose the center's fractional parts to be 0.5.")
    print(f"Minimal cells = 8 * R - 2 * (2 + 2) = 8 * {R} - 8")
    print(f"Minimal cells = {8*R} - 8 = {min_cells}")
    
    print(f"\nThe minimal and maximal numbers of cells are ({min_cells}, {max_cells}).")

solve_grid_crossing()