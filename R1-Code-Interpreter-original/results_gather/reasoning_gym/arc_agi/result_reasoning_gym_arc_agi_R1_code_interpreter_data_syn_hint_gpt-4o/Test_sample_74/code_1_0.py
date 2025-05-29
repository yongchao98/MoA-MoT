import numpy as np

def find_largest_number_section(grid):
    # Convert the grid to a numpy array for easier manipulation
    grid = np.array(grid)
    
    # Find the largest number in the grid
    largest_number = np.max(grid)
    
    # Find all positions of the largest number
    positions = np.argwhere(grid == largest_number)
    
    # Determine the bounding box for these positions
    min_row, min_col = np.min(positions, axis=0)
    max_row, max_col = np.max(positions, axis=0)
    
    # Extract the section of the grid defined by the bounding box
    section = grid[min_row:max_row+1, min_col:max_col+1]
    
    return section

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 9, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 9, 4, 9, 9, 9, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 9, 9, 9, 9, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 9, 9, 4, 9, 9, 9, 9, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 9, 4, 4, 9, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 6, 4, 4, 4, 6, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 6, 4, 4, 6, 4, 6, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 6, 4, 6, 4, 6, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4, 4],
    [4, 4, 4, 6, 4, 6, 4, 6, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 5, 4, 4],
    [4, 4, 6, 4, 6, 6, 4, 6, 4, 4, 4, 4, 4, 4, 4, 5, 5, 4, 5, 5, 4, 4],
    [4, 4, 6, 4, 4, 4, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
]

# Get the output grid
output_grid = find_largest_number_section(input_grid)

# Print the output grid
print(output_grid)