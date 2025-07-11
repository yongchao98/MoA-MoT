import numpy as np

def create_initial_grid(pattern_type):
    """
    Creates the initial 8x8 grid with marked cells based on the specified pattern.
    - pattern_type=1: (i+j) even -> Top-Left, odd -> Bottom-Right
    - pattern_type=2: (i+j) even -> Bottom-Right, odd -> Top-Left
    """
    grid = np.zeros((8, 8), dtype=int)
    # Iterate through the 4x4 grid of 2x2 sub-grids
    for i in range(4):
        for j in range(4):
            # Top-left coordinate of the 2x2 sub-grid in the 8x8 grid
            subgrid_row, subgrid_col = 2 * i, 2 * j
            is_checkerboard_even = (i + j) % 2 == 0

            if (pattern_type == 1 and is_checkerboard_even) or (pattern_type == 2 and not is_checkerboard_even):
                # Mark the top-left corner of the 2x2 sub-grid
                grid[subgrid_row, subgrid_col] = 1
            else:
                # Mark the bottom-right corner of the 2x2 sub-grid
                grid[subgrid_row + 1, subgrid_col + 1] = 1
    return grid

def apply_transformations(grid):
    """
    Applies the specified transformations:
    1. Reflection over the line y=x (transpose)
    2. Rotation 90 degrees clockwise
    """
    # 1. Reflection over y=x is a transpose operation
    reflected_grid = np.transpose(grid)

    # 2. Rotation 90 degrees clockwise
    # This can be achieved by transposing and then reversing each row.
    # Note: The problem asks for rotation of the *reflected* grid.
    rotated_grid = np.fliplr(reflected_grid.T)
    return rotated_grid

def count_subgrids(grid, subgrid_size, target_sum):
    """
    Counts the number of sub-grids of a given size with a specific sum of marked cells.
    """
    count = 0
    grid_size = grid.shape[0]
    num_subgrids_per_dim = grid_size - subgrid_size + 1
    
    for r in range(num_subgrids_per_dim):
        for c in range(num_subgrids_per_dim):
            subgrid = grid[r : r + subgrid_size, c : c + subgrid_size]
            if np.sum(subgrid) == target_sum:
                count += 1
    return count

# --- Main execution ---
# Analyze the first checkerboard pattern
initial_grid_1 = create_initial_grid(1)
final_grid_1 = apply_transformations(initial_grid_1)
count_pattern_1 = count_subgrids(final_grid_1, 4, 2)

# Analyze the second checkerboard pattern
initial_grid_2 = create_initial_grid(2)
final_grid_2 = apply_transformations(initial_grid_2)
count_pattern_2 = count_subgrids(final_grid_2, 4, 2)

# The result is the maximum of the two counts
max_count = max(count_pattern_1, count_pattern_2)

print(f"Number of 4x4 sub-grids with exactly two marked cells for pattern 1: {count_pattern_1}")
print(f"Number of 4x4 sub-grids with exactly two marked cells for pattern 2: {count_pattern_2}")
print(f"The maximum number of 4x4 sub-grids with exactly two marked cells is: {max_count}")

<<<9>>>