import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation and counting problem.
    """
    
    max_count = 0
    
    # There are two possible interpretations of the checkerboard pattern.
    # We test both and find the maximum result.
    for case in ['A', 'B']:
        # 1. & 2. Construct the Initial Grid
        grid = np.zeros((8, 8), dtype=int)
        
        # Populate the grid based on the checkerboard pattern
        for sr in range(4): # sub-grid row index
            for sc in range(4): # sub-grid column index
                # Top-left corner of the 2x2 sub-grid
                base_r, base_c = 2 * sr, 2 * sc
                
                is_even = (sr + sc) % 2 == 0
                
                if (case == 'A' and is_even) or (case == 'B' and not is_even):
                    # Mark top-left corner
                    grid[base_r, base_c] = 1
                else:
                    # Mark bottom-right corner
                    grid[base_r + 1, base_c + 1] = 1

        # 3. Apply Transformations
        # Reflect over the line y=x (transpose)
        reflected_grid = grid.T
        
        # Rotate 90 degrees clockwise
        # In numpy, k=-1 specifies a clockwise rotation
        final_grid = np.rot90(reflected_grid, k=-1)
        
        # 4. Count Sub-grids with Exactly Two Marked Cells
        current_case_count = 0
        subgrid_size = 4
        
        for r in range(grid.shape[0] - subgrid_size + 1):
            for c in range(grid.shape[1] - subgrid_size + 1):
                # Extract the 4x4 sub-grid
                subgrid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
                
                # Check if it has exactly two marked cells
                if np.sum(subgrid) == 2:
                    current_case_count += 1
        
        if current_case_count > max_count:
            max_count = current_case_count

    print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is: {max_count}")

solve_grid_problem()