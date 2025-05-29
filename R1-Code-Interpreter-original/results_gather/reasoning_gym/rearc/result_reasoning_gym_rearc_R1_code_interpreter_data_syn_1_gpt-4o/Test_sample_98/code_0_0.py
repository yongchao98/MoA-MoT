def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        output_row = [8 if x == 3 else x for x in row]
        output_grid.append(output_row)
    return output_grid

input_grid = [
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 4, 1, 1],
    [1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 3, 1, 1, 1],
    [1, 1, 1, 3, 2, 1, 1]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))