import numpy as np

def solve_grid_problem():
    """
    This function solves the grid problem by creating the grid,
    applying transformations, and counting the sub-grids that meet the criteria.
    """
    # Step 1: Construct the Initial Grid
    N = 8
    initial_grid = np.zeros((N, N), dtype=int)

    # The 8x8 grid is composed of a 4x4 arrangement of 2x2 sub-grids.
    for i in range(N // 2):  # Sub-grid row index (0 to 3)
        for j in range(N // 2):  # Sub-grid column index (0 to 3)
            # Use checkerboard pattern for the sub-grids
            if (i + j) % 2 == 0:
                # Mark top-left corner of the 2x2 sub-grid
                r, c = 2 * i, 2 * j
                initial_grid[r, c] = 1
            else:
                # Mark bottom-right corner of the 2x2 sub-grid
                r, c = 2 * i + 1, 2 * j + 1
                initial_grid[r, c] = 1

    # Step 2: Apply Transformations
    # Reflect over y=x (transpose)
    reflected_grid = initial_grid.transpose()

    # Rotate 90 degrees clockwise
    final_grid = np.rot90(reflected_grid, k=-1)

    # Step 3 & 4: Count sub-grids with exactly two marked cells
    subgrid_size = 4
    num_subgrids_per_dim = N - subgrid_size + 1
    
    indicators = []
    
    for r in range(num_subgrids_per_dim):
        for c in range(num_subgrids_per_dim):
            sub_grid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            marked_cells = np.sum(sub_grid)
            
            # Check if the sub-grid has exactly two marked cells
            if marked_cells == 2:
                indicators.append(1)
            else:
                indicators.append(0)

    # Step 5: Final Calculation and Output
    total_count = sum(indicators)
    
    print("The final calculation is the sum of indicators for each of the 25 possible 4x4 sub-grids.")
    print("Each indicator is 1 if the sub-grid has exactly two marked cells, and 0 otherwise.")
    
    # Format the equation string as requested
    equation_str = " + ".join(map(str, indicators))
    print(f"\n{equation_str} = {total_count}")
    
    print(f"\nThe maximum number of 4x4 sub-grids with exactly two marked cells is {total_count}.")

solve_grid_problem()
<<<0>>>