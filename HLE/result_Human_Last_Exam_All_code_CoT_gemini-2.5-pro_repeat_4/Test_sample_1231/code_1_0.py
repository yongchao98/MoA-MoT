import numpy as np

def solve():
    """
    Solves the grid transformation and counting problem.
    """
    # Step 1 & 2: Create the initial 8x8 grid and place marked cells.
    initial_grid = np.zeros((8, 8), dtype=int)
    for i in range(4):  # Index for 2x2 sub-grids
        for j in range(4):  # Index for 2x2 sub-grids
            if (i + j) % 2 == 0:
                # Mark top-left corner of the 2x2 sub-grid
                r, c = 2 * i, 2 * j
                initial_grid[r, c] = 1
            else:
                # Mark bottom-right corner of the 2x2 sub-grid
                r, c = 2 * i + 1, 2 * j + 1
                initial_grid[r, c] = 1

    # Step 3: Apply the combined transformation to get the final grid.
    # A point (r, c) after reflection over y=x becomes (c, r).
    # Then, rotating (c, r) 90 degrees clockwise results in (r, 7-c).
    final_grid = np.zeros((8, 8), dtype=int)
    for r in range(8):
        for c in range(8):
            if initial_grid[r, c] == 1:
                new_r, new_c = r, 7 - c
                final_grid[new_r, new_c] = 1

    # Step 4: Iterate through all 4x4 sub-grids and count those with exactly two marks.
    subgrids_with_two_marks = 0
    subgrid_size = 4
    grid_size = 8
    
    # There are (8 - 4 + 1) x (8 - 4 + 1) = 25 possible sub-grids.
    for r_start in range(grid_size - subgrid_size + 1):
        for c_start in range(grid_size - subgrid_size + 1):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r_start:r_start + subgrid_size, c_start:c_start + subgrid_size]
            
            # Count the number of marked cells in the sub-grid
            marked_cells_count = np.sum(sub_grid)
            
            # Check if the count is exactly 2
            if marked_cells_count == 2:
                subgrids_with_two_marks += 1
    
    # Step 5: Print the final result.
    print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is:")
    print(subgrids_with_two_marks)

solve()