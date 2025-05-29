def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        # Duplicate each element in the row
        new_row = []
        for element in row:
            new_row.extend([element, element])
        output_grid.append(new_row)
    return output_grid

# Test input grid
input_grid = [
    [5, 5],
    [5, 5]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))