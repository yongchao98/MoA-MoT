def transform_grid(input_grid):
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    output_grid = input_grid.copy()
    
    while len(output_grid) < num_cols:
        output_grid.extend(input_grid)
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3],
    [3, 3]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))