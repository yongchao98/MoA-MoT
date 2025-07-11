import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation and counting problem.
    """
    size = 8
    
    # Step 1 & 2: Create the initial 8x8 grid with marked cells
    initial_grid = np.zeros((size, size), dtype=int)
    num_subgrids_axis = size // 2

    for i in range(num_subgrids_axis):
        for j in range(num_subgrids_axis):
            if (i + j) % 2 == 0:
                # Top-left corner of the 2x2 sub-grid
                initial_grid[2 * i, 2 * j] = 1
            else:
                # Bottom-right corner of the 2x2 sub-grid
                initial_grid[2 * i + 1, 2 * j + 1] = 1

    # Step 3: Apply transformations
    # Transformation 1: Reflect over y=x (transpose)
    reflected_grid = np.transpose(initial_grid)

    # Transformation 2: Rotate 90 degrees clockwise
    # np.rot90(m, k=1) is a 90-degree counter-clockwise rotation.
    # k=3 gives a 270-degree counter-clockwise rotation, which is a 90-degree clockwise rotation.
    final_grid = np.rot90(reflected_grid, k=3)

    # Step 4 & 5: Count sub-grids with exactly two marked cells
    subgrid_size = 4
    count_two_marked = 0
    
    # There are (8 - 4 + 1) * (8 - 4 + 1) = 25 possible 4x4 sub-grids
    num_subgrids_total = (size - subgrid_size + 1)**2
    
    for r in range(size - subgrid_size + 1):
        for c in range(size - subgrid_size + 1):
            subgrid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            if np.sum(subgrid) == 2:
                count_two_marked += 1
                
    # Step 6: Print the result.
    # The problem asks for "the maximum number", which refers to counting these subgrids
    # in the single resulting configuration.
    # The "equation" is simply the final count.
    print(f"The total number of 4x4 sub-grids is 25.")
    print(f"After checking all 25 sub-grids in the final configuration, we found that the number of sub-grids with exactly two marked cells is:")
    print(count_two_marked)

solve_grid_problem()
<<<0>>>