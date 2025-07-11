def solve_grid_crossing():
    """
    Calculates the minimal and maximal number of grid cells a circle can cross.
    """
    # Radius of the circle
    R = 500

    # The number of cells a simple closed curve crosses (C) is the sum of
    # its horizontal (H) and vertical (V) grid line crossings.
    # C = H + V

    # Number of horizontal lines the circle's y-span covers.
    # The interval is (y0 - R, y0 + R), with length 2*R.
    # Since R is an integer and tangency is forbidden, the number of integers
    # in this interval is always 2*R.
    num_horizontal_lines = 2 * R

    # Each horizontal line is crossed twice.
    H = 2 * num_horizontal_lines

    # The same logic applies to vertical lines.
    num_vertical_lines = 2 * R

    # Each vertical line is crossed twice.
    V = 2 * num_vertical_lines

    # The total number of cells crossed is C = H + V = 8 * R.
    total_cells = H + V

    # Since the result is independent of the circle's specific valid position,
    # the minimum and maximum values are the same.
    min_cells = total_cells
    max_cells = total_cells
    
    # --- Outputting the calculation steps ---
    print(f"The radius of the circle is R = {R}.")
    print(f"The number of horizontal grid lines crossed is N_h = 2 * {R} = {num_horizontal_lines}.")
    print(f"The total number of horizontal crossings is H = 2 * {num_horizontal_lines} = {H}.")
    print(f"The number of vertical grid lines crossed is N_v = 2 * {R} = {num_vertical_lines}.")
    print(f"The total number of vertical crossings is V = 2 * {num_vertical_lines} = {V}.")
    print(f"The total number of cells crossed is C = H + V = {H} + {V} = {total_cells}.")
    print(f"\nSince the number of crossed cells is constant, the minimal and maximal numbers are the same.")
    
    # --- Final Answer ---
    print(f"The final answer in the form (min, max) is ({min_cells}, {max_cells}).")

solve_grid_crossing()