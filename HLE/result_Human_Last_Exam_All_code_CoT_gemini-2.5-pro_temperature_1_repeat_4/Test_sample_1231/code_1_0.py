import numpy as np

def create_initial_grid(config_type):
    """
    Creates the initial 8x8 grid based on the configuration type.
    config_type=1: top-left for even i+j, bottom-right for odd i+j.
    config_type=2: bottom-right for even i+j, top-left for odd i+j.
    """
    grid = np.zeros((8, 8), dtype=int)
    for i in range(4):
        for j in range(4):
            is_even = (i + j) % 2 == 0
            if (config_type == 1 and is_even) or (config_type == 2 and not is_even):
                # Mark top-left corner of the 2x2 sub-grid
                grid[2 * i, 2 * j] = 1
            else:
                # Mark bottom-right corner of the 2x2 sub-grid
                grid[2 * i + 1, 2 * j + 1] = 1
    return grid

def apply_transformations(grid):
    """
    Reflects the grid over y=x and then rotates 90 degrees clockwise.
    This is equivalent to reversing each row.
    """
    # Create an empty grid for the final state
    final_grid = np.zeros((8, 8), dtype=int)
    
    # A point (r, c) after reflection over y=x is (c, r).
    # A point (r, c) after 90-deg clockwise rotation is (c, 7-r).
    # Combined: (r, c) -> (c, r) -> (r, 7-c).
    # So, final_grid[r][7-c] = grid[r][c]
    # which is equivalent to reversing each row.
    for r in range(8):
        final_grid[r, :] = grid[r, ::-1]
        
    return final_grid

def count_subgrids_with_two_marks(grid):
    """
    Counts the number of 4x4 sub-grids with exactly two marked cells.
    """
    count = 0
    # The top-left corner of a 4x4 sub-grid can be from (0,0) to (4,4)
    for r_start in range(5):
        for c_start in range(5):
            # Extract the 4x4 sub-grid
            sub_grid = grid[r_start:r_start+4, c_start:c_start+4]
            # Sum the marked cells
            if np.sum(sub_grid) == 2:
                count += 1
    return count

def solve_puzzle():
    """
    Solves the puzzle by checking both configurations and finding the maximum.
    """
    # Configuration 1
    initial_grid_1 = create_initial_grid(1)
    final_grid_1 = apply_transformations(initial_grid_1)
    count1 = count_subgrids_with_two_marks(final_grid_1)
    
    # Configuration 2
    initial_grid_2 = create_initial_grid(2)
    final_grid_2 = apply_transformations(initial_grid_2)
    count2 = count_subgrids_with_two_marks(final_grid_2)
    
    max_count = max(count1, count2)
    
    print(f"Number of 4x4 sub-grids with exactly two marked cells for Configuration 1: {count1}")
    print(f"Number of 4x4 sub-grids with exactly two marked cells for Configuration 2: {count2}")
    print(f"The maximum number is max({count1}, {count2}) = {max_count}")

# Run the solver
solve_puzzle()