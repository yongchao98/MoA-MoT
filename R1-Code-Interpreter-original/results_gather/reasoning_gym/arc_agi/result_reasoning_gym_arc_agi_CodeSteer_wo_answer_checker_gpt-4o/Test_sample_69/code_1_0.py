def extract_output_grid(input_grid):
    # Determine the size of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # The output grid is typically a 3x3 section
    output_size = 3
    
    # Identify the starting point for the output grid
    # Based on the examples, it seems to be located in the lower part of the grid
    start_row = rows - output_size - 2
    start_col = cols - output_size - 2
    
    # Extract the output grid
    output_grid = []
    for i in range(start_row, start_row + output_size):
        output_grid.append(input_grid[i][start_col:start_col + output_size])
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3],
    [3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3],
    [3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3]
]

# Extract the output grid
output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))