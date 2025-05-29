def transform_grid(input_grid):
    # Find the rows with non-zero numbers
    non_zero_rows = [row for row in input_grid if any(cell != 0 for cell in row)]
    
    # Calculate how many times we need to repeat the non-zero rows to fill the grid
    repeat_count = len(input_grid) // len(non_zero_rows)
    
    # Construct the output grid by repeating the non-zero rows
    output_grid = []
    for _ in range(repeat_count):
        output_grid.extend(non_zero_rows)
    
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