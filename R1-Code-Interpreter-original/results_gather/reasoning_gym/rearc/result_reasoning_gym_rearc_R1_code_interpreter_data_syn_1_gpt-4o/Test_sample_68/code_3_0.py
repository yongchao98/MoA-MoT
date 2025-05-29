import numpy as np

def extract_grid(input_grid):
    # Convert the input grid to a numpy array for easier manipulation
    grid = np.array(input_grid)
    
    # Determine the border values
    border_values = set(grid[0]) | set(grid[-1]) | set(grid[:, 0]) | set(grid[:, -1])
    
    # Flatten the grid and count occurrences of each value, excluding border values
    flat_grid = grid.flatten()
    value_counts = {}
    for value in flat_grid:
        if value not in border_values:
            if value in value_counts:
                value_counts[value] += 1
            else:
                value_counts[value] = 1
    
    # Check if value_counts is empty and handle it
    if not value_counts:
        raise ValueError("No non-border values found.")
    
    # Find the most frequent non-border value
    prominent_value = max(value_counts, key=value_counts.get)
    
    # Find the region of the prominent value
    rows, cols = np.where(grid == prominent_value)
    min_row, max_row = min(rows), max(rows)
    min_col, max_col = min(cols), max(cols)
    
    # Extract the grid around the prominent value
    # Adjust the size based on the examples, here we assume a 10x6 grid for consistency
    extracted_grid = grid[min_row:max_row+1, min_col:max_col+1]
    
    return extracted_grid

# Test input grid
input_grid = [
    [8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
]

# Extract the grid
output_grid = extract_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))