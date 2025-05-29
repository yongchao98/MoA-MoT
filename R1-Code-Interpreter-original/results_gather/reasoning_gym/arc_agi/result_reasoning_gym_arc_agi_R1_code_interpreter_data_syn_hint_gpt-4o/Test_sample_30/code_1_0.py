# Define the test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 0, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 1, 9, 9, 9, 0, 0, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 1, 9, 1, 9, 9, 0, 9, 9, 9],
    [9, 9, 9, 8, 9, 8, 9, 9, 1, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 8, 9, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 5, 6, 6, 6, 8, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 5, 9, 5, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 5, 5, 5, 9, 9, 9, 9],
    [6, 7, 6, 6, 6, 0, 6, 9, 9, 9, 5, 9, 9, 9, 9, 9],
    [6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9]
]

# Function to extract a 7x7 subgrid based on specific patterns
def extract_subgrid(grid):
    # Identify key rows and columns based on specific values
    key_values = [8, 0, 1]
    start_row, start_col = None, None
    
    # Find the first occurrence of key values in the grid
    for i, row in enumerate(grid):
        for j, value in enumerate(row):
            if value in key_values:
                start_row, start_col = i, j
                break
        if start_row is not None:
            break
    
    # Ensure the subgrid is within bounds
    if start_row is not None and start_col is not None:
        start_row = max(0, min(start_row, len(grid) - 7))
        start_col = max(0, min(start_col, len(grid[0]) - 7))
        subgrid = [row[start_col:start_col+7] for row in grid[start_row:start_row+7]]
        return subgrid
    else:
        return []

# Extract the subgrid
output_grid = extract_subgrid(input_grid)

# Print the output grid
print(output_grid)