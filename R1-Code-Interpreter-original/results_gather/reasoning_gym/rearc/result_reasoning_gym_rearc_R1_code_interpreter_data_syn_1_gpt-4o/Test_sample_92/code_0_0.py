def transform_grid(input_grid):
    # Copy the input grid to output grid
    output_grid = [row[:] for row in input_grid]
    
    # Transform the last column to 1
    for row in output_grid:
        row[-1] = 1
    
    # Transform the third column to 2
    for row in output_grid:
        row[2] = 2
    
    # Transform the first column to 3
    for row in output_grid:
        row[0] = 3
    
    return output_grid

# Test input grid
input_grid = [
    [5, 8, 5, 5, 5],
    [5, 8, 5, 5, 5],
    [5, 8, 5, 5, 9],
    [5, 8, 7, 5, 9],
    [1, 8, 7, 5, 9]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))