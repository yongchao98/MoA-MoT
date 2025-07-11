import numpy as np

def calculate_subgrids_with_two_marks(initial_grid):
    """
    Transforms the grid and counts the number of 4x4 sub-grids with exactly two marked cells.
    
    The transformation is a reflection over y=x followed by a 90-degree clockwise rotation.
    This combined operation is equivalent to a vertical flip of the grid.
    """
    
    # Vertically flip the grid to perform the combined transformation
    final_grid = np.flipud(initial_grid)
    
    count_of_two = 0
    # There are (8-4+1) x (8-4+1) = 25 possible 4x4 sub-grids
    # Their top-left corners can be at any row/col from 0 to 4
    for r_start in range(5):
        for c_start in range(5):
            # Extract the 4x4 sub-grid
            sub_grid = final_grid[r_start:r_start+4, c_start:c_start+4]
            
            # Sum the marked cells in the sub-grid
            marked_cells = np.sum(sub_grid)
            
            # Check if the count is exactly two
            if marked_cells == 2:
                count_of_two += 1
                
    return count_of_two

def main():
    """
    Generates the two possible initial grid patterns, calculates the result for each,
    and prints the maximum.
    """
    
    # Pattern A: (R+C) even -> top-left, (R+C) odd -> bottom-right
    grid_A = np.zeros((8, 8), dtype=int)
    # Pattern B: (R+C) even -> bottom-right, (R+C) odd -> top-left
    grid_B = np.zeros((8, 8), dtype=int)

    # Iterate through the 4x4 array of 2x2 sub-grids
    for R in range(4):  # Sub-grid row index
        for C in range(4):  # Sub-grid column index
            # Get the top-left coordinate (r, c) of the 2x2 sub-grid
            r, c = 2 * R, 2 * C
            
            if (R + C) % 2 == 0:  # 'Even' sub-grids on the checkerboard
                grid_A[r, c] = 1              # Pattern A: Mark top-left
                grid_B[r + 1, c + 1] = 1      # Pattern B: Mark bottom-right
            else:  # 'Odd' sub-grids on the checkerboard
                grid_A[r + 1, c + 1] = 1      # Pattern A: Mark bottom-right
                grid_B[r, c] = 1              # Pattern B: Mark top-left

    # Calculate the number of sub-grids with two marks for each pattern
    result_A = calculate_subgrids_with_two_marks(grid_A)
    result_B = calculate_subgrids_with_two_marks(grid_B)
    
    # The answer is the maximum of the two results
    max_subgrids = max(result_A, result_B)
    
    print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is: {max_subgrids}")

if __name__ == "__main__":
    main()