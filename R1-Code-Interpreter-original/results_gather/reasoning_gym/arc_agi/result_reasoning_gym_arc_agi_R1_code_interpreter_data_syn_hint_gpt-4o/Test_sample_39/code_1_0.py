def transform_grid(input_grid):
    output_grid = []
    for i, row in enumerate(input_grid):
        new_row = row[:]
        # Determine which column to replace with 5 based on the row index
        new_row[i % 3] = 5
        output_grid.append(new_row)
    return output_grid

input_grid = [
    [0, 9, 9],
    [9, 0, 9],
    [9, 9, 0],
    [0, 9, 9],
    [9, 9, 0],
    [9, 0, 9],
    [0, 9, 9],
    [9, 9, 0],
    [9, 0, 9],
    [9, 9, 0],
    [0, 9, 9],
    [9, 0, 9],
    [0, 9, 9],
    [9, 9, 0],
    [9, 0, 9],
    [0, 9, 9],
    [9, 0, 9]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))