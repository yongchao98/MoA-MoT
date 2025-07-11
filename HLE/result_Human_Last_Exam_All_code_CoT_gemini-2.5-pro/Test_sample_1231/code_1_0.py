import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation and counting problem.
    """
    N = 8

    # 1. & 2. Create and populate the initial grid
    initial_grid = np.zeros((N, N), dtype=int)
    for i in range(N // 2):
        for j in range(N // 2):
            # Use checkerboard pattern for 2x2 sub-grids
            if (i + j) % 2 == 0:
                # Mark top-left corner of the 2x2 sub-grid
                r, c = 2 * i, 2 * j
                initial_grid[r, c] = 1
            else:
                # Mark bottom-right corner of the 2x2 sub-grid
                r, c = 2 * i + 1, 2 * j + 1
                initial_grid[r, c] = 1

    # 3. Apply transformations to get the final grid
    final_grid = np.zeros((N, N), dtype=int)
    for r in range(N):
        for c in range(N):
            if initial_grid[r, c] == 1:
                # Reflection over y=x: (r, c) -> (c, r)
                # 90-deg clockwise rotation: (r', c') -> (c', N-1-r')
                # Combined: (r, c) -> (r, N-1-c)
                new_r, new_c = r, N - 1 - c
                final_grid[new_r, new_c] = 1

    # 4. & 5. & 6. Count 4x4 sub-grids with exactly two marked cells
    subgrid_size = 4
    num_subgrids_with_two_marks = 0
    
    # There are (N - subgrid_size + 1) x (N - subgrid_size + 1) possible sub-grids
    for r_start in range(N - subgrid_size + 1):
        for c_start in range(N - subgrid_size + 1):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r_start : r_start + subgrid_size, c_start : c_start + subgrid_size]
            
            # Count marked cells in the sub-grid
            marked_cells_count = np.sum(sub_grid)
            
            # Check if the count is exactly 2
            if marked_cells_count == 2:
                num_subgrids_with_two_marks += 1
    
    # Print the final result
    final_count = num_subgrids_with_two_marks
    print(f"The number of 4x4 sub-grids with exactly two marked cells is: {final_count}")

solve_grid_problem()