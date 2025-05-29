def transform_grid(input_grid):
    output_grid = []
    for row in input_grid:
        new_row = row[:]
        # Find the last occurrence of zero and replace it with 5
        for j in range(len(row) - 1, -1, -1):
            if row[j] == 0:
                new_row[j] = 5
                break
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

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))