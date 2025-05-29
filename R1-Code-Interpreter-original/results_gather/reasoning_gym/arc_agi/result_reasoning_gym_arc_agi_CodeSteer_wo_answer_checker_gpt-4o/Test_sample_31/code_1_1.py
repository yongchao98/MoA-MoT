def extract_output_grid(input_grid):
    # Determine the dimensions of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Assume the border value is the value at the top-left corner
    border_value = input_grid[0][0]
    
    # Initialize variables to find the bounding box of non-border values
    min_row, max_row = rows, 0
    min_col, max_col = cols, 0
    
    # Iterate over the grid to find the bounds of the non-border values
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != border_value:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    # Extract the subgrid based on the identified bounds
    output_grid = []
    for r in range(min_row, max_row + 1):
        output_grid.append(input_grid[r][min_col:max_col + 1])
    
    return output_grid

# Test input grid
input_grid = [
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 9, 9, 9, 9, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 2, 2, 9, 9, 9, 9, 2, 2, 2, 2, 6, 6],
    [6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 3, 3, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
]

# Get the output grid
output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))