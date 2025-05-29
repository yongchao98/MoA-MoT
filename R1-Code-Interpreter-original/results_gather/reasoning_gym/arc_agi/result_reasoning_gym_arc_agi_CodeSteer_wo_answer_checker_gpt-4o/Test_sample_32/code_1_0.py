def transform_grid(input_grid):
    # Remove columns with '5'
    output_grid = []
    for row in input_grid:
        new_row = [num for i, num in enumerate(row) if i % 4 != 3]
        output_grid.append(new_row)
    return output_grid

# Test input grid
input_grid = [
    [6, 6, 4, 5, 6, 6, 4, 5, 6, 6, 4, 5, 6, 6, 4],
    [6, 4, 4, 5, 6, 4, 4, 5, 6, 4, 4, 5, 6, 4, 4],
    [9, 9, 9, 5, 4, 4, 4, 5, 4, 4, 4, 5, 4, 4, 4],
    [6, 9, 9, 5, 6, 9, 9, 5, 6, 4, 4, 5, 6, 4, 4],
    [6, 6, 9, 5, 6, 6, 9, 5, 6, 6, 9, 5, 6, 6, 4],
    [6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9],
    [9, 9, 9, 5, 9, 9, 9, 5, 9, 9, 9, 5, 9, 9, 9],
    [6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9],
    [6, 6, 9, 5, 6, 6, 9, 5, 6, 6, 9, 5, 6, 6, 9],
    [6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9],
    [9, 9, 9, 5, 9, 9, 9, 5, 9, 9, 9, 5, 9, 9, 9],
    [6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9, 5, 6, 9, 9],
    [6, 6, 9, 5, 6, 6, 9, 5, 6, 6, 9, 5, 6, 6, 9]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))