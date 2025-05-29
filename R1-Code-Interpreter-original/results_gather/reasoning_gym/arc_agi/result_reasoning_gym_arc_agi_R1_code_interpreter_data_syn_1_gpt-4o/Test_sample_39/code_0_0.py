def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        new_row = []
        for i in range(len(row)):
            if i > 0 and i < len(row) - 1 and row[i-1] == 9 and row[i] == 0 and row[i+1] == 9:
                new_row.append(5)
            else:
                new_row.append(row[i])
        output_grid.append(new_row)
    return output_grid

# Test input grid
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