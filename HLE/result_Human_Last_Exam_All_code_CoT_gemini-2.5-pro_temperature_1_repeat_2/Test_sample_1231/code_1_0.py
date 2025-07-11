import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation problem by simulating the steps.
    """
    # Step 1: Initialize an 8x8 grid
    # Using numpy for easier array manipulation, initialized with zeros.
    initial_grid = np.zeros((8, 8), dtype=int)

    # Step 2: Populate the initial grid based on the checkerboard pattern
    # The 8x8 grid is composed of a 4x4 arrangement of 2x2 sub-grids.
    # Let (i, j) be the coordinates of the 2x2 sub-grid (from 0 to 3).
    print("Generating the initial 8x8 grid with marked cells...")
    for i in range(4):
        for j in range(4):
            # The choice of corner depends on a checkerboard pattern of the sub-grids.
            if (i + j) % 2 == 0:
                # If (i+j) is even, mark the top-left corner of the 2x2 sub-grid.
                # The top-left corner's coordinates are (2*i, 2*j).
                initial_grid[2 * i, 2 * j] = 1
            else:
                # If (i+j) is odd, mark the bottom-right corner.
                # The bottom-right corner's coordinates are (2*i+1, 2*j+1).
                initial_grid[2 * i + 1, 2 * j + 1] = 1
    
    # print("Initial Grid:")
    # print(initial_grid)

    # Step 3 & 4: Apply transformations to get the final grid
    # A cell at (r, c) is reflected over y=x to (c, r).
    # Then, (c, r) is rotated 90-deg clockwise to (r, 8-1-c) = (r, 7-c).
    final_grid = np.zeros((8, 8), dtype=int)
    print("Applying transformations (reflection then rotation)...")
    for r in range(8):
        for c in range(8):
            if initial_grid[r, c] == 1:
                new_r, new_c = r, 7 - c
                final_grid[new_r, new_c] = 1

    # print("\nFinal Grid after transformations:")
    # print(final_grid)

    # Step 5 & 6: Count 4x4 sub-grids with exactly two marked cells
    # An 8x8 grid has (8-4+1)x(8-4+1) = 5x5 = 25 possible 4x4 sub-grids.
    num_subgrids_with_two_marks = 0
    print("\nCounting marked cells in all 25 possible 4x4 sub-grids...")
    
    # Iterate through all possible top-left corners (r_start, c_start) of a 4x4 sub-grid.
    for r_start in range(5):
        for c_start in range(5):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r_start:r_start + 4, c_start:c_start + 4]
            
            # Count the number of marked cells (1s) in the sub-grid
            marked_cells_count = np.sum(sub_grid)
            
            # Check if the count is exactly two
            if marked_cells_count == 2:
                num_subgrids_with_two_marks += 1

    print(f"\nFinal calculation:")
    print(f"The number of 4x4 sub-grids with exactly two marked cells is: {num_subgrids_with_two_marks}")

solve_grid_problem()
<<<0>>>