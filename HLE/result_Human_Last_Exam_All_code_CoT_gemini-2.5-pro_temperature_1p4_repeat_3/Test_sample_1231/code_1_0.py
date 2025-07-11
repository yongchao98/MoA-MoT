import numpy as np

def solve_grid_problem():
    """
    Solves the described grid transformation and counting problem.
    """
    # 1. Construct the Initial Grid
    grid = np.zeros((8, 8), dtype=int)
    for i in range(4):
        for j in range(4):
            # The choice of which corner to mark depends on a checkerboard pattern
            # on the 4x4 grid of 2x2 sub-grids.
            if (i + j) % 2 == 0:
                # Mark top-left corner of the 2x2 sub-grid
                grid[2 * i, 2 * j] = 1
            else:
                # Mark bottom-right corner of the 2x2 sub-grid
                grid[2 * i + 1, 2 * j + 1] = 1

    # 2. Apply Transformations
    # a) Reflection over y=x is a transpose.
    reflected_grid = np.transpose(grid)

    # b) Rotate 90 degrees clockwise. np.rot90 rotates counter-clockwise,
    # so we use k=-1 for one clockwise rotation.
    final_grid = np.rot90(reflected_grid, k=-1)

    # 3. Count 4x4 sub-grids with exactly two marked cells
    count = 0
    grid_size = 8
    subgrid_size = 4
    
    # There are (8-4+1) x (8-4+1) = 25 possible sub-grids
    for r in range(grid_size - subgrid_size + 1):
        for c in range(grid_size - subgrid_size + 1):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            
            # Sum the marked cells in the sub-grid
            marked_cells = np.sum(sub_grid)
            
            if marked_cells == 2:
                count += 1
    
    # 4. Print the final answer
    # The problem asks for "each number in the final equation"
    # The equation is simply the total count.
    print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is:")
    print(count)

solve_grid_problem()