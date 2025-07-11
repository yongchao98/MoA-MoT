import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation problem and finds the number of 4x4 sub-grids
    with exactly two marked cells.
    """
    # 1. & 2. Construct the Initial Grid and Mark the Cells
    # Initialize an 8x8 grid with zeros
    grid = np.zeros((8, 8), dtype=int)

    # Iterate through the 2x2 sub-grids' top-left corners
    # The sub-grids form a 4x4 meta-grid, indexed by (i, j)
    for i in range(4):
        for j in range(4):
            # Use a checkerboard pattern for marking
            if (i + j) % 2 == 0:
                # Mark the top-left cell of the 2x2 sub-grid
                grid[2 * i, 2 * j] = 1
            else:
                # Mark the bottom-right cell of the 2x2 sub-grid
                grid[2 * i + 1, 2 * j + 1] = 1

    # 3. Apply Transformations
    # Reflection over the line y=x is a transpose operation
    reflected_grid = grid.transpose()

    # Rotation 90 degrees clockwise
    # np.rot90 with k=-1 rotates 90 degrees clockwise
    final_grid = np.rot90(reflected_grid, k=-1)
    
    # 4. Count Sub-grids with exactly two marked cells
    subgrid_count = 0
    # There are (8-4+1) x (8-4+1) = 25 possible 4x4 sub-grids
    grid_size = 8
    subgrid_size = 4
    
    for r in range(grid_size - subgrid_size + 1):
        for c in range(grid_size - subgrid_size + 1):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            
            # Count the marked cells in the sub-grid
            marked_cells = np.sum(sub_grid)
            
            # Check if the count is exactly 2
            if marked_cells == 2:
                subgrid_count += 1
                
    # 5. Find the Result
    print(f"The initial grid is:\n{grid}\n")
    print(f"The final grid after reflection and rotation is:\n{final_grid}\n")
    print(f"The number of 4x4 sub-grids with exactly two marked cells is: {subgrid_count}")

solve_grid_problem()
<<<0>>>