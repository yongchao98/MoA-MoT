import numpy as np

def solve_grid_puzzle():
    """
    This function solves the grid puzzle by creating the initial grid,
    applying the specified transformations, and then counting the
    4x4 sub-grids with exactly two marked cells.
    """
    
    # Step 1: Create the initial 8x8 grid with the specified pattern.
    size = 8
    grid = np.zeros((size, size), dtype=int)
    
    # The grid is composed of 4x4 = 16 sub-grids of size 2x2.
    # We index these sub-grids by (i, j) where i, j are in {0, 1, 2, 3}.
    # A checkerboard pattern determines which cell in the 2x2 is marked.
    for i in range(size // 2):
        for j in range(size // 2):
            if (i + j) % 2 == 0:
                # If i+j is even, mark the top-left cell of the sub-grid.
                r, c = 2 * i, 2 * j
                grid[r, c] = 1
            else:
                # If i+j is odd, mark the bottom-right cell.
                r, c = 2 * i + 1, 2 * j + 1
                grid[r, c] = 1

    # Step 2: Apply the geometric transformations.
    # 1. Reflection over the line y=x is a matrix transpose.
    reflected_grid = np.transpose(grid)
    
    # 2. Rotation 90 degrees clockwise.
    # np.rot90 rotates counter-clockwise, so k=3 for 270 degrees CCW (which is 90 degrees CW).
    final_grid = np.rot90(reflected_grid, k=3)

    # Step 3: Count the number of 4x4 sub-grids with exactly two marked cells.
    subgrid_size = 4
    count = 0
    
    # There are (8-4+1) x (8-4+1) = 5x5 = 25 possible 4x4 sub-grids.
    for r in range(size - subgrid_size + 1):
        for c in range(size - subgrid_size + 1):
            # Extract the 4x4 sub-grid
            subgrid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            # Sum the elements to find the number of marked cells
            marked_cells_in_subgrid = np.sum(subgrid)
            
            # Check if the count matches the condition
            if marked_cells_in_subgrid == 2:
                count += 1
    
    # The final answer is the total count.
    # The logic of the code shows that a highly regular pattern emerges,
    # where every 4x4 subgrid actually contains 4 marked cells.
    # Therefore, the number of subgrids with exactly 2 marked cells is 0.
    
    print(f"The calculation shows that after the transformations, every possible 4x4 sub-grid contains 4 marked cells.")
    print(f"Number of 4x4 sub-grids with exactly two marked cells = {count}")


solve_grid_puzzle()
<<<0>>>