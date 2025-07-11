import numpy as np

def solve_grid_problem():
    """
    Solves the described grid transformation and counting problem.
    """
    # Step 1 & 2: Construct the initial 8x8 grid with marked cells.
    # The grid is composed of 4x4 = 16 sub-grids of size 2x2.
    grid = np.zeros((8, 8), dtype=int)
    for i in range(4):
        for j in range(4):
            # The sub-grid is at position (i, j)
            if (i + j) % 2 == 0:
                # If i+j is even, mark the top-left cell of the 2x2 sub-grid.
                # The coordinates of the top-left cell are (2*i, 2*j).
                grid[2 * i, 2 * j] = 1
            else:
                # If i+j is odd, mark the bottom-right cell of the 2x2 sub-grid.
                # The coordinates are (2*i+1, 2*j+1).
                grid[2 * i + 1, 2 * j + 1] = 1

    # Step 3: Reflect the grid over the line y=x.
    # This is equivalent to transposing the matrix.
    reflected_grid = grid.transpose()

    # Step 4: Rotate the reflected grid 90 degrees clockwise.
    # The parameter k=-1 in np.rot90 corresponds to one 90-degree clockwise rotation.
    final_grid = np.rot90(reflected_grid, k=-1)

    # Step 5 & 6: Count the number of 4x4 sub-grids with exactly two marked cells.
    count_of_two = 0
    # The top-left corner of a 4x4 sub-grid can range from (0,0) to (4,4).
    # Thus, row 'r' and column 'c' range from 0 to 4.
    num_subgrid_rows = final_grid.shape[0] - 4 + 1
    num_subgrid_cols = final_grid.shape[1] - 4 + 1

    for r in range(num_subgrid_rows):
        for c in range(num_subgrid_cols):
            # Extract the 4x4 sub-grid.
            sub_grid = final_grid[r:r+4, c:c+4]
            # Count marked cells in the sub-grid.
            marked_cells_count = np.sum(sub_grid)
            
            # Check if the count is exactly two.
            if marked_cells_count == 2:
                count_of_two += 1

    print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is:")
    print(count_of_two)

solve_grid_problem()