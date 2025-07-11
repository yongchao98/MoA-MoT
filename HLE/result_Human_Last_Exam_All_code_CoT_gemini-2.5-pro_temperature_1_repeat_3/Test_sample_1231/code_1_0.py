import numpy as np

def solve_grid_problem():
    """
    This function solves the grid transformation and counting problem.
    """
    # Step 1: Construct the initial 8x8 grid.
    # The grid is composed of 2x2 sub-grids.
    # The sub-grid at index (i, j) where 0 <= i, j <= 3
    # has a marked cell at the top-left if i+j is even,
    # and at the bottom-right if i+j is odd.
    initial_grid = np.zeros((8, 8), dtype=int)
    for i in range(4):
        for j in range(4):
            if (i + j) % 2 == 0:
                # Top-left corner of the (i,j)-th 2x2 sub-grid
                row, col = 2 * i, 2 * j
                initial_grid[row, col] = 1
            else:
                # Bottom-right corner of the (i,j)-th 2x2 sub-grid
                row, col = 2 * i + 1, 2 * j + 1
                initial_grid[row, col] = 1

    # Step 2: Apply the transformations.
    # First, reflect over the line y=x, which is a matrix transpose.
    reflected_grid = np.transpose(initial_grid)

    # Second, rotate the result 90 degrees clockwise.
    # In numpy, rot90 with k=-1 rotates clockwise.
    final_grid = np.rot90(reflected_grid, k=-1)

    # Step 3 & 4: Iterate through all 4x4 sub-grids and check the sum.
    # An 8x8 grid has (8-4+1) x (8-4+1) = 25 possible 4x4 sub-grids.
    # The top-left corner of a sub-grid can be at (r, c) where 0 <= r, c <= 4.
    count_of_two_marked_cells = 0
    for r in range(5):
        for c in range(5):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r:r+4, c:c+4]
            
            # Count the number of marked cells (sum of elements)
            marked_cells_sum = np.sum(sub_grid)
            
            # Check if the sum is exactly 2
            if marked_cells_sum == 2:
                count_of_two_marked_cells += 1

    # Step 5: Print the final result.
    print(f"The initial grid has {np.sum(initial_grid)} marked cells.")
    print(f"The final grid also has {np.sum(final_grid)} marked cells.")
    print(f"Number of 4x4 sub-grids checked: 25")
    print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is: {count_of_two_marked_cells}")

solve_grid_problem()