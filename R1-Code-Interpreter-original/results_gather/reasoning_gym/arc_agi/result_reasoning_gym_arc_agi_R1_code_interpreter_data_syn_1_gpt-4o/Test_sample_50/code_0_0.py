def predict_output(input_grid):
    output_grid = []
    for i in range(0, 9, 3):
        row = []
        for j in range(0, 9, 3):
            subgrid = [input_grid[x][j:j+3] for x in range(i, i+3)]
            if any(5 in subrow for subrow in subgrid):
                row.append(6)
            else:
                row.append(3)
        output_grid.append(row)
    return output_grid

input_grid = [
    [3, 3, 3, 2, 3, 5, 3, 2, 3, 3, 3],
    [5, 3, 3, 2, 3, 3, 3, 2, 3, 5, 3],
    [3, 3, 5, 2, 3, 3, 3, 2, 3, 3, 3],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [3, 3, 3, 2, 5, 3, 3, 2, 3, 3, 5],
    [3, 3, 5, 2, 3, 5, 3, 2, 3, 3, 3],
    [3, 3, 3, 2, 3, 3, 3, 2, 3, 5, 3],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3],
    [3, 5, 3, 2, 3, 3, 5, 2, 3, 3, 3],
    [3, 3, 3, 2, 3, 3, 3, 2, 5, 3, 5]
]

output_grid = predict_output(input_grid)
print(output_grid)