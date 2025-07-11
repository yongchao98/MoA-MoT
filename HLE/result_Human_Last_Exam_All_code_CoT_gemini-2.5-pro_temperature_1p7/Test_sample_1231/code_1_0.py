import numpy as np

def solve_grid_problem():
    """
    This function solves the grid problem by modeling the setup,
    applying the transformations, and counting the resulting sub-grids.
    """

    # The problem implies we can choose the starting pattern of the checkerboard.
    # We will test both possibilities and take the maximum result.
    max_count = 0

    # There are two possible checkerboard patterns.
    for pattern_choice in [0, 1]:
        # 1. Construct the Initial Grid
        initial_grid = np.zeros((8, 8), dtype=int)
        print(f"--- Testing Checkerboard Pattern {pattern_choice + 1} ---")
        print("Step 1: Construct the initial 8x8 grid.")
        for i in range(4): # 2x2 sub-grid row index
            for j in range(4): # 2x2 sub-grid column index
                if (i + j) % 2 == pattern_choice:
                    # Mark top-left cell of the 2x2 sub-grid
                    r, c = 2 * i, 2 * j
                    initial_grid[r, c] = 1
                    print(f"Sub-grid({i},{j}) has i+j={i+j}. Marking top-left at ({r},{c})")
                else:
                    # Mark bottom-right cell of the 2x2 sub-grid
                    r, c = 2 * i + 1, 2 * j + 1
                    initial_grid[r, c] = 1
                    print(f"Sub-grid({i},{j}) has i+j={i+j}. Marking bottom-right at ({r},{c})")

        # 2. Apply Transformations
        # 2a. Reflect over the line y=x (transpose)
        print("\nStep 2: Apply transformations.")
        reflected_grid = initial_grid.T
        print("Reflected grid over y=x (transposed).")

        # 2b. Rotate 90 degrees clockwise
        final_grid = np.rot90(reflected_grid, k=-1)
        print("Rotated the result 90 degrees clockwise.")
        print("Final Grid Configuration:")
        print(final_grid)

        # 3. Count 4x4 sub-grids with exactly two marked cells
        print("\nStep 3: Count 4x4 sub-grids with exactly two marked cells.")
        count_with_two_marked = 0
        subgrid_size = 4
        grid_size = 8
        
        # There are (8-4+1) x (8-4+1) = 25 possible sub-grids
        for r_start in range(grid_size - subgrid_size + 1):
            for c_start in range(grid_size - subgrid_size + 1):
                sub_grid = final_grid[r_start:r_start+subgrid_size, c_start:c_start+subgrid_size]
                marked_cells_count = np.sum(sub_grid)
                if marked_cells_count == 2:
                    count_with_two_marked += 1
        
        print(f"\nFor this pattern, number of 4x4 sub-grids with exactly two marked cells is: {count_with_two_marked}")
        
        if count_with_two_marked > max_count:
            max_count = count_with_two_marked

    print("\n--- Final Conclusion ---")
    print(f"The maximum number of 4x4 sub-grids that can be formed with exactly two marked cells is the maximum count found across all possible patterns.")
    print(f"Final Answer = {max_count}")

solve_grid_problem()