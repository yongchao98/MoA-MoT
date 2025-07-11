import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation puzzle.
    """
    # Step 1 & 2: Create the initial 8x8 grid and place the marks.
    initial_grid = np.zeros((8, 8), dtype=int)
    # Iterate over the 4x4 grid of 2x2 sub-grids
    for i in range(4):
        for j in range(4):
            # The top-left corner of the 2x2 sub-grid in the 8x8 grid
            row, col = i * 2, j * 2
            # Use a checkerboard pattern for placing marks
            if (i + j) % 2 == 0:
                # Mark top-left corner of the sub-grid
                initial_grid[row, col] = 1
            else:
                # Mark bottom-right corner of the sub-grid
                initial_grid[row + 1, col + 1] = 1

    # Step 3: Apply transformations
    # Reflection over the line y=x is a transpose
    reflected_grid = initial_grid.T
    
    # Rotation 90 degrees clockwise. In numpy, k=-1 does a clockwise rotation.
    final_grid = np.rot90(reflected_grid, k=-1)

    # Step 4: Count 4x4 sub-grids with exactly two marked cells
    count = 0
    # A 4x4 sub-grid can start at any row from 0 to 4 and any col from 0 to 4
    num_rows, num_cols = final_grid.shape
    subgrid_size = 4
    for r in range(num_rows - subgrid_size + 1):
        for c in range(num_cols - subgrid_size + 1):
            sub_grid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            # The sum of the sub-grid equals the number of marked cells
            if np.sum(sub_grid) == 2:
                count += 1
    
    # Step 5: Output the result
    print(f"{count}")

solve_grid_problem()
<<<0>>>