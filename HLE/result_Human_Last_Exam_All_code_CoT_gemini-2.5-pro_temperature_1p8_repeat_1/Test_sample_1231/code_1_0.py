import numpy as np

def solve_grid_problem():
    """
    Solves the grid problem by creating the grid, applying transformations,
    and counting the 4x4 sub-grids with exactly two marked cells.
    """
    # Step 1: Construct the Initial Grid
    # An 8x8 grid composed of 2x2 sub-grids.
    grid_size = 8
    initial_grid = np.zeros((grid_size, grid_size), dtype=int)

    # Populate the grid based on the checkerboard pattern of marked cells.
    # A 2x2 sub-grid's top-left corner is at (2*I, 2*J).
    for i in range(grid_size):
        for j in range(grid_size):
            # Determine the index (I, J) of the 2x2 sub-grid.
            I, J = i // 2, j // 2
            
            # The pattern alternates based on the sum of sub-grid indices.
            if (I + J) % 2 == 0:
                # For even sums, mark the top-left cell of the 2x2 sub-grid.
                if i % 2 == 0 and j % 2 == 0:
                    initial_grid[i, j] = 1
            else:
                # For odd sums, mark the bottom-right cell.
                if i % 2 == 1 and j % 2 == 1:
                    initial_grid[i, j] = 1

    # Step 2: Apply Geometric Transformations
    # Reflection over y=x (transpose) followed by a 90-degree clockwise
    # rotation is equivalent to a horizontal flip.
    transformed_grid = np.fliplr(initial_grid)

    # Step 3 & 4: Count Marked Cells in Sub-grids
    subgrid_size = 4
    count_of_two_marked = 0
    
    # The number of possible sub-grids along one dimension.
    num_subgrids_dim = grid_size - subgrid_size + 1

    # Iterate through all possible top-left corners (r, c) of a 4x4 sub-grid.
    for r in range(num_subgrids_dim):
        for c in range(num_subgrids_dim):
            # Extract the 4x4 sub-grid.
            sub_grid = transformed_grid[r:r+subgrid_size, c:c+subgrid_size]
            
            # Sum the marked cells in the sub-grid.
            if np.sum(sub_grid) == 2:
                count_of_two_marked += 1

    # Step 5: Print the Final Answer
    print("The size of the main grid is {0}x{0}.".format(grid_size))
    print("The size of the sub-grids to check is {0}x{0}.".format(subgrid_size))
    print(f"The total number of possible {subgrid_size}x{subgrid_size} sub-grids is {num_subgrids_dim * num_subgrids_dim}.")
    print("After all transformations, we find the number of sub-grids with exactly 2 marked cells.")
    print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is: {count_of_two_marked}")

solve_grid_problem()