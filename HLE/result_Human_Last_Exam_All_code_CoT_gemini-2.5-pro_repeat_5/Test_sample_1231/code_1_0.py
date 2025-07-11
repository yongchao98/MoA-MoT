import collections

def solve_grid_problem():
    """
    Solves the described grid transformation and counting problem.
    """
    # Step 1: Construct the initial 8x8 grid with marked cells
    initial_grid = [[0] * 8 for _ in range(8)]
    initial_marked_cells = []
    for i in range(4):
        for j in range(4):
            # The 2x2 sub-grid's top-left corner is at (2*i, 2*j)
            # The choice of marked cell depends on the checkerboard pattern of the sub-grids
            if (i + j) % 2 == 0:
                # Mark the top-left cell of the sub-grid
                r, c = 2 * i, 2 * j
            else:
                # Mark the bottom-right cell of the sub-grid
                r, c = 2 * i + 1, 2 * j + 1
            initial_grid[r][c] = 1
            initial_marked_cells.append((r, c))

    # Step 2: Apply the transformations to find the final grid configuration
    final_grid = [[0] * 8 for _ in range(8)]
    final_marked_cells = []
    for r, c in initial_marked_cells:
        # Transformation 1: Reflection over y=x -> (r, c) becomes (c, r)
        # Transformation 2: Rotate 90 degrees clockwise -> (c, r) becomes (r, 7-c)
        final_r, final_c = r, 7 - c
        final_grid[final_r][final_c] = 1
        final_marked_cells.append((final_r, final_c))

    # Step 3 & 4: Count the number of 4x4 sub-grids with exactly two marked cells
    subgrids_with_two_marks = 0
    grid_size = 8
    subgrid_size = 4
    
    # There are (grid_size - subgrid_size + 1) possible start positions for rows and columns
    for start_r in range(grid_size - subgrid_size + 1):
        for start_c in range(grid_size - subgrid_size + 1):
            marked_count = 0
            # Iterate through the cells of the 4x4 sub-grid
            for r_offset in range(subgrid_size):
                for c_offset in range(subgrid_size):
                    r = start_r + r_offset
                    c = start_c + c_offset
                    if final_grid[r][c] == 1:
                        marked_count += 1
            
            # Check if the count is exactly two
            if marked_count == 2:
                subgrids_with_two_marks += 1
    
    # Output the result. The problem asks for the final number, which is the count.
    # The wording "final equation" is interpreted as providing the final calculated value.
    print(f"The number of 4x4 sub-grids with exactly two marked cells is: {subgrids_with_two_marks}")

solve_grid_problem()