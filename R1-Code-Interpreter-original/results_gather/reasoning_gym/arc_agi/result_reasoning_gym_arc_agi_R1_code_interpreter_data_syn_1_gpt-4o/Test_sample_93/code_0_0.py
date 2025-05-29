def repeat_grid_vertically(input_grid):
    return input_grid + input_grid

input_grid = [
    [7, 6, 6, 6, 6],
    [6, 7, 1, 6, 6],
    [6, 1, 9, 9, 6],
    [6, 1, 6, 8, 8]
]

output_grid = repeat_grid_vertically(input_grid)
for row in output_grid:
    print(row)