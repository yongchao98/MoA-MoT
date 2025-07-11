import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation problem and finds the number of 4x4 sub-grids
    with exactly two marked cells.
    """
    # Step 1: Construct the initial 8x8 grid.
    initial_grid = np.zeros((8, 8), dtype=int)
    for r in range(8):
        for c in range(8):
            # Determine the 2x2 sub-grid indices (i, j)
            i = r // 2
            j = c // 2
            # Apply the checkerboard pattern rule
            if (i + j) % 2 == 0:  # Even sub-grid: mark top-left
                if r == 2 * i and c == 2 * j:
                    initial_grid[r, c] = 1
            else:  # Odd sub-grid: mark bottom-right
                if r == 2 * i + 1 and c == 2 * j + 1:
                    initial_grid[r, c] = 1

    # Step 2: Apply the geometric transformations.
    # a) Reflect over the line y=x (transpose)
    reflected_grid = np.transpose(initial_grid)
    
    # b) Rotate 90 degrees clockwise
    # In NumPy, rot90 with k=-1 rotates 90 degrees clockwise.
    final_grid = np.rot90(reflected_grid, k=-1)

    # Step 3: Count the 4x4 sub-grids with exactly two marked cells.
    count = 0
    # There are (8-4+1) * (8-4+1) = 25 possible 4x4 sub-grids.
    # The top-left corner of a sub-grid can be from (0,0) to (4,4).
    for r_start in range(5):
        for c_start in range(5):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r_start:r_start+4, c_start:c_start+4]
            # Sum the marked cells in the sub-grid
            marked_cells = np.sum(sub_grid)
            # Check if the count is exactly 2
            if marked_cells == 2:
                count += 1
    
    # Step 4: Output the result.
    print(count)

solve_grid_problem()