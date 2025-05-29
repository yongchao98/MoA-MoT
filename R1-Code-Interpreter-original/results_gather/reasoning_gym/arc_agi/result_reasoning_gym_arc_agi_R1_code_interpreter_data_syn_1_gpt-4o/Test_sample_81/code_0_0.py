def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        sixes = [x for x in row if x == 6]
        non_sixes = [x for x in row if x != 6]
        new_row = sixes + non_sixes
        output_grid.append(new_row)
    return output_grid

input_grid = [
    [6, 2, 6, 2, 2],
    [8, 6, 6, 8, 6],
    [6, 6, 9, 6, 6],
    [1, 6, 6, 1, 6],
    [3, 6, 6, 6, 6]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(row)