def extend_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    # Extend numbers vertically downwards
    for col in range(cols):
        for row in range(rows):
            if input_grid[row][col] != 0:
                num = input_grid[row][col]
                # Extend downwards
                for r in range(row + 1, rows):
                    if input_grid[r][col] == 0:
                        output_grid[r][col] = num
                    else:
                        break

    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 5, 9, 9, 9, 9, 9, 9, 9, 9, 2],
    [9, 9, 9, 9, 2, 5, 9, 9, 9, 9, 9, 9, 9],
    [9, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 5, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 1, 1, 1],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 1, 1, 1],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 1, 1, 1],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]
]

output_grid = extend_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))