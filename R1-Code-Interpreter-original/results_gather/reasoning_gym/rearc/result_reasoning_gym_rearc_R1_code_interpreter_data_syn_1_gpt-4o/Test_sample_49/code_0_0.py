def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        if any(row):  # If there's any non-zero number in the row
            output_grid.append([0] * len(row))
        else:  # If the row is entirely zero
            output_grid.append([5] * len(row))
    return output_grid

# Test input grid
input_grid = [
    [0, 9, 0],
    [0, 0, 0],
    [0, 9, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))