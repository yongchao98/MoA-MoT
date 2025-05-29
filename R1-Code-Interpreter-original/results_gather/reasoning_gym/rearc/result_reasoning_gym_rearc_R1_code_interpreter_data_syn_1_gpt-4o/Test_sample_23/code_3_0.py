def transform_grid(input_grid):
    # Find the rows with non-zero numbers
    non_zero_rows = [row for row in input_grid if any(cell != 0 for cell in row)]
    
    # Construct the output grid by repeating the non-zero rows
    output_grid = []
    while len(output_grid) < len(input_grid):
        output_grid.extend(non_zero_rows)
    
    # Trim the output grid to match the input grid's number of rows
    output_grid = output_grid[:len(input_grid)]
    
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