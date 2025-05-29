def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        output_row = [8 if x == 9 else 0 for x in row]
        output_grid.append(output_row)
    return output_grid

# Test input grid
input_grid = [
    [0, 3, 6, 0, 9],
    [3, 0, 6, 9, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))