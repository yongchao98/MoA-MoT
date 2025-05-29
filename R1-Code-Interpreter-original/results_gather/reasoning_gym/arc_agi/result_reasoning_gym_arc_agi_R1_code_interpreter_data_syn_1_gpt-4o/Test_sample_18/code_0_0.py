def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        output_row = []
        for value in row:
            if value == 7:
                output_row.append(7)
            else:
                output_row.append(9)
        output_grid.append(output_row)
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 2, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 2, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [7, 7, 7, 9, 9, 9, 9, 9, 9],
    [7, 9, 7, 9, 9, 9, 9, 9, 9],
    [9, 7, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))