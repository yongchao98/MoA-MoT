def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        new_row = [1 if x == 5 else x for x in row]
        output_grid.append(new_row)
    return output_grid

input_grid = [
    [1, 1, 8, 5],
    [6, 1, 0, 5],
    [6, 0, 2, 4],
    [9, 8, 1, 1]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))