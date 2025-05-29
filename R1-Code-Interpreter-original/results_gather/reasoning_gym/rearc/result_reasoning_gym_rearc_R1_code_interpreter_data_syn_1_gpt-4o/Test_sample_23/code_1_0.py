def transform_grid(input_grid):
    # Find the first row with non-zero numbers
    non_zero_row = None
    for row in input_grid:
        if any(cell != 0 for cell in row):
            non_zero_row = row
            break
    
    # If a non-zero row is found, repeat it to fill the grid
    if non_zero_row:
        output_grid = [non_zero_row] * len(input_grid)
    else:
        output_grid = input_grid  # If no non-zero row is found, return the input grid
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 6],
    [7, 7, 7, 7, 7, 6],
    [7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))