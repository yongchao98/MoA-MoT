def extract_subgrid(input_grid):
    # Determine the dimensions of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the non-background value range
    min_row, max_row = rows, 0
    min_col, max_col = cols, 0
    
    # Background value is assumed to be the most common value in the border
    background_value = input_grid[0][0]
    
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != background_value:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    # Extract the subgrid
    subgrid = []
    for r in range(min_row, max_row + 1):
        subgrid.append(input_grid[r][min_col:max_col + 1])
    
    return subgrid

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

# Extract and print the output grid
output_grid = extract_subgrid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))