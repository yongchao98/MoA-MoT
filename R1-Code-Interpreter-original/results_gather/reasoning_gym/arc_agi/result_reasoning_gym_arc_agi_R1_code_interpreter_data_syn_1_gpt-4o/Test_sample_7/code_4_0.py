def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        if row.count(9) > row.count(4):
            output_grid.append([9] * len(row))
        else:
            output_row = []
            for i in range(len(row)):
                if i % 3 == 0:
                    output_row.append(9)
                else:
                    output_row.append(4)
            output_grid.append(output_row)
    return output_grid

input_grid = [
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, 4, 4]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))