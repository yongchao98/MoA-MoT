def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to the output grid

    # Shift the sequence '9 9 9' to the right by 3 positions
    for i in range(len(input_grid)):
        if input_grid[i][3:6] == [9, 9, 9]:
            output_grid[i][3:6] = [0, 0, 0]
            output_grid[i][6:9] = [9, 9, 9]

    # Shift the sequence '9 6 9' to the right by 3 positions
    for i in range(len(input_grid)):
        if input_grid[i][3:6] == [9, 6, 9]:
            output_grid[i][3:6] = [0, 0, 0]
            output_grid[i][6:9] = [9, 6, 9]

    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 9, 9, 9, 0],
    [0, 0, 0, 9, 6, 9, 0],
    [0, 0, 0, 9, 9, 9, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 6, 0, 0]
]

output_grid = transform_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))