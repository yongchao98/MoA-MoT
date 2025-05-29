def transform_grid(input_grid):
    output_grid = []
    for i, row in enumerate(input_grid):
        new_row = row[:]
        if i % 3 == 0:
            new_row[0] = 5
        elif i % 3 == 1:
            new_row[1] = 5
        else:
            new_row[2] = 5
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