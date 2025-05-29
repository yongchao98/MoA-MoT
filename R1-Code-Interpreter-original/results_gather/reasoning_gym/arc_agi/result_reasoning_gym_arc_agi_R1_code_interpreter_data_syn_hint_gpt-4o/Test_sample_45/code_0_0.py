def extract_subgrid(input_grid):
    # Find the bounds of the non-zero subgrid
    min_row, max_row = len(input_grid), 0
    min_col, max_col = len(input_grid[0]), 0
    
    for i, row in enumerate(input_grid):
        for j, value in enumerate(row):
            if value != 0:
                min_row = min(min_row, i)
                max_row = max(max_row, i)
                min_col = min(min_col, j)
                max_col = max(max_col, j)
    
    # Extract the subgrid
    output_grid = []
    for i in range(min_row, max_row + 1):
        output_grid.append(input_grid[i][min_col:max_col + 1])
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 5, 2, 0, 0],
    [0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]

# Get the output grid
output_grid = extract_subgrid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))