def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        if all(x == 9 for x in row):
            output_grid.append([9] * len(row))
        else:
            transformed_row = []
            for i in range(len(row)):
                if i % 3 == 0:
                    transformed_row.append(9)
                else:
                    transformed_row.append(4)
            output_grid.append(transformed_row)
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