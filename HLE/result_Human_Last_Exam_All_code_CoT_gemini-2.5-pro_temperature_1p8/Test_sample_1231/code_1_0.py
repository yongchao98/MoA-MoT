import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation and counting problem.
    """
    # Step 1 & 2: Construct the initial 8x8 grid and place marked cells.
    # The grid is tiled by 2x2 sub-grids, indexed by (i, j) where i,j in [0,3].
    initial_grid = np.zeros((8, 8), dtype=int)
    for i in range(4):
        for j in range(4):
            # A checkerboard pattern on the sub-grid indices (i, j)
            # determines the position of the marked cell.
            if (i + j) % 2 == 0:
                # If (i+j) is even, mark the top-left cell of the sub-grid.
                row, col = 2 * i, 2 * j
                initial_grid[row, col] = 1
            else:
                # If (i+j) is odd, mark the bottom-right cell.
                row, col = 2 * i + 1, 2 * j + 1
                initial_grid[row, col] = 1

    # Step 3: Apply the geometric transformations.

    # First, reflect the grid over the line y=x. This is equivalent to
    # transposing the matrix.
    reflected_grid = np.transpose(initial_grid)

    # Second, rotate the reflected grid 90 degrees clockwise.
    # The parameter k=-1 corresponds to a -90 degree (clockwise) rotation.
    final_grid = np.rot90(reflected_grid, k=-1)

    # Step 4: Iterate through all possible 4x4 sub-grids and count those
    # with exactly two marked cells.
    # An 8x8 grid has (8-4+1) x (8-4+1) = 25 possible 4x4 sub-grids.
    num_subgrids_with_two_marks = 0
    grid_size = 8
    subgrid_size = 4
    
    # The starting top-left corner (r, c) of a subgrid can range from
    # (0, 0) to (4, 4).
    for r in range(grid_size - subgrid_size + 1):
        for c in range(grid_size - subgrid_size + 1):
            # Extract the 4x4 sub-grid from the final grid.
            sub_grid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            
            # Count the number of marked cells (1s) in the sub-grid.
            marked_cells_count = np.sum(sub_grid)
            
            # Check if the count is exactly two.
            if marked_cells_count == 2:
                num_subgrids_with_two_marks += 1
    
    # Step 5: Print the final result.
    print(f"{num_subgrids_with_two_marks}")

solve_grid_problem()