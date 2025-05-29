def transform_grid(input_grid):
    # Remove rows that are entirely zeros
    return [row for row in input_grid if any(cell != 0 for cell in row)]

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 9, 0, 0, 4, 0, 4, 0, 0, 5, 0, 0, 6, 6, 6],
    [0, 0, 9, 0, 4, 0, 4, 0, 5, 5, 5, 0, 0, 6, 0],
    [9, 9, 0, 0, 4, 0, 4, 0, 0, 5, 0, 0, 6, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 9, 9, 0, 4, 4, 0, 0, 5, 0, 5, 0, 0, 6, 0],
    [0, 9, 0, 0, 0, 0, 4, 0, 5, 0, 5, 0, 6, 6, 6],
    [9, 9, 0, 0, 4, 4, 0, 0, 5, 0, 5, 0, 0, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))