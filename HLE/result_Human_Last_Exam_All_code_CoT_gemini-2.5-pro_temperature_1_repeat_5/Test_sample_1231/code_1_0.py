import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation and counting problem.
    """
    N = 8
    
    # 1. Construct the Initial Grid
    # The grid is 8x8, composed of 4x4 = 16 sub-grids of size 2x2.
    # We use 1 for a marked cell and 0 for an unmarked cell.
    initial_grid = np.zeros((N, N), dtype=int)
    
    # Iterate over the 2x2 sub-grids
    for i in range(N // 2):  # sub-grid row index
        for j in range(N // 2):  # sub-grid column index
            # Top-left corner of the 2x2 sub-grid
            top_left_r, top_left_c = 2 * i, 2 * j
            
            # Apply checkerboard pattern to choose which cell to mark
            if (i + j) % 2 == 0:
                # Mark the top-left cell of the 2x2 sub-grid
                initial_grid[top_left_r, top_left_c] = 1
            else:
                # Mark the bottom-right cell of the 2x2 sub-grid
                initial_grid[top_left_r + 1, top_left_c + 1] = 1

    # 2. Apply the Transformations
    # The transformation consists of a reflection over y=x, then a 90-degree clockwise rotation.
    # Let's trace a point (r, c):
    # - Reflection over y=x: (r, c) -> (c, r)
    # - 90-degree clockwise rotation around the center: (c, r) -> (r, N-1-c)
    # So, the final position of the cell at initial_grid[r][c] is transformed_grid[r][N-1-c].
    # This is equivalent to a horizontal flip of the grid.
    transformed_grid = np.zeros((N, N), dtype=int)
    for r in range(N):
        for c in range(N):
            transformed_grid[r, c] = initial_grid[r, N - 1 - c]

    # 3. Count marked cells in all 4x4 sub-grids
    # An 8x8 grid has (8-4+1)x(8-4+1) = 5x5 = 25 possible 4x4 sub-grids.
    count_of_two_marked_cells = 0
    subgrid_size = 4
    
    # Iterate over all possible top-left corners of 4x4 sub-grids
    for r_start in range(N - subgrid_size + 1):
        for c_start in range(N - subgrid_size + 1):
            # Extract the 4x4 sub-grid
            sub_grid = transformed_grid[r_start:r_start + subgrid_size, c_start:c_start + subgrid_size]
            
            # Calculate the number of marked cells in the sub-grid
            marked_cells_sum = np.sum(sub_grid)
            
            # Check if the sum is exactly 2
            if marked_cells_sum == 2:
                count_of_two_marked_cells += 1
                
    # 4. Final Answer
    # The problem asks for the maximum number. Since the initial configuration is fixed
    # (up to a symmetric choice that yields the same result), the answer is the calculated count.
    print("The maximum number of 4x4 sub-grids with exactly two marked cells is:")
    print(count_of_two_marked_cells)

solve_grid_problem()