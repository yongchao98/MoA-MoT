import numpy as np

def solve_grid_problem():
    """
    Solves the grid transformation and sub-grid counting problem.
    The function considers two possible initial checkerboard patterns,
    applies the specified transformations, and counts the number of
    4x4 sub-grids with exactly two marked cells in each case.
    The final answer is the maximum of these two counts.
    """

    results = {}

    # There are two possible starting checkerboard patterns.
    # Pattern A: (0,0) sub-grid has its top-left cell marked.
    # Pattern B: (0,0) sub-grid has its bottom-right cell marked.
    for pattern_type in ['A', 'B']:
        
        # 1. Construct the Initial 8x8 Grid
        initial_grid = np.zeros((8, 8), dtype=int)
        # Iterate over the 4x4 array of 2x2 sub-grids
        for i in range(4):  # sub-grid row index from 0 to 3
            for j in range(4):  # sub-grid column index from 0 to 3
                # Top-left corner of the sub-grid (i, j)
                subgrid_row, subgrid_col = 2 * i, 2 * j
                
                # Determine mark position based on checkerboard pattern
                is_even = (i + j) % 2 == 0
                
                # Apply Pattern A or B
                if (pattern_type == 'A' and is_even) or \
                   (pattern_type == 'B' and not is_even):
                    # Mark the top-left cell of the sub-grid
                    initial_grid[subgrid_row, subgrid_col] = 1
                else:
                    # Mark the bottom-right cell of the sub-grid
                    initial_grid[subgrid_row + 1, subgrid_col + 1] = 1

        # 2. Apply Transformations
        # Step 2a: Reflect over the line y=x. This is a matrix transpose.
        # Both initial grid patterns are symmetric across the y=x diagonal,
        # so the transpose does not change the grid.
        reflected_grid = initial_grid.transpose()

        # Step 2b: Rotate 90 degrees clockwise.
        # np.rot90(m, k=-1) rotates 90 degrees clockwise.
        final_grid = np.rot90(reflected_grid, k=-1)

        # 3. Count 4x4 Sub-grids with Exactly Two Marked Cells
        count_with_2_marks = 0
        # Iterate through all possible top-left corners of a 4x4 sub-grid
        for r in range(5):  # r from 0 to 4
            for c in range(5):  # c from 0 to 4
                # Extract the 4x4 sub-grid
                subgrid = final_grid[r:r + 4, c:c + 4]
                # Sum the elements to find the number of marked cells
                if np.sum(subgrid) == 2:
                    count_with_2_marks += 1
        
        results[pattern_type] = count_with_2_marks

    # 4. Determine the Maximum count
    max_count = max(results.values())
    
    print(f"For Pattern A (top-left first), the number of 4x4 sub-grids with 2 marks is: {results['A']}")
    print(f"For Pattern B (bottom-right first), the number of 4x4 sub-grids with 2 marks is: {results['B']}")
    print(f"The maximum number of such sub-grids is {max_count}.")

solve_grid_problem()