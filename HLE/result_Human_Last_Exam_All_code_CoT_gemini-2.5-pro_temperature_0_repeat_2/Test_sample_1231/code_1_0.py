def solve_grid_problem():
    """
    This function solves the problem by following these steps:
    1. Constructs the initial 8x8 grid with the specified pattern of marked cells.
    2. Applies the reflection and rotation transformations.
    3. Iterates through all 4x4 sub-grids of the final grid.
    4. Counts how many of these sub-grids have exactly two marked cells.
    5. Prints the final count.
    """
    grid_size = 8
    
    # Step 1: Construct the Initial Grid
    initial_grid = [[0] * grid_size for _ in range(grid_size)]
    
    num_subgrids_dim = grid_size // 2
    for i in range(num_subgrids_dim):
        for j in range(num_subgrids_dim):
            row_start, col_start = 2 * i, 2 * j
            
            # Determine which cell to mark based on the checkerboard pattern
            # of the 2x2 sub-grids.
            if (i + j) % 2 == 0:
                # Mark top-left corner for "white" squares of the sub-grid checkerboard
                initial_grid[row_start][col_start] = 1
            else:
                # Mark bottom-right corner for "black" squares
                initial_grid[row_start + 1][col_start + 1] = 1

    # Step 2: Apply Transformations
    # a) Reflect over the line y=x (transpose)
    reflected_grid = [[0] * grid_size for _ in range(grid_size)]
    for r in range(grid_size):
        for c in range(grid_size):
            reflected_grid[c][r] = initial_grid[r][c]

    # b) Rotate the reflected grid 90 degrees clockwise
    final_grid = [[0] * grid_size for _ in range(grid_size)]
    for r in range(grid_size):
        for c in range(grid_size):
            final_grid[c][grid_size - 1 - r] = reflected_grid[r][c]

    # Step 3 & 4: Count 4x4 sub-grids with exactly two marked cells
    count_with_two_marked = 0
    subgrid_size = 4
    
    # Iterate through all possible top-left corners of 4x4 sub-grids
    for r_start in range(grid_size - subgrid_size + 1):
        for c_start in range(grid_size - subgrid_size + 1):
            
            marked_cells_in_subgrid = 0
            # Sum the marked cells in the current 4x4 sub-grid
            for r in range(r_start, r_start + subgrid_size):
                for c in range(c_start, c_start + subgrid_size):
                    if final_grid[r][c] == 1:
                        marked_cells_in_subgrid += 1
            
            # If the sum is 2, we've found a sub-grid that meets the condition
            if marked_cells_in_subgrid == 2:
                count_with_two_marked += 1

    # The problem asks for the maximum number. Since the grid and transformations
    # are fixed, there is only one possible outcome for the number of such sub-grids.
    print(f"{count_with_two_marked}")

solve_grid_problem()
<<<0>>>