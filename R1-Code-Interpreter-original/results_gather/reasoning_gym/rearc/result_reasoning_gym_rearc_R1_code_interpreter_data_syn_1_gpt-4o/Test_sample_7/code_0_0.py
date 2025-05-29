def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy input grid to output grid

    # Identify the number that forms a vertical or horizontal line
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 0:
                num = input_grid[i][j]
                break
        else:
            continue
        break

    # Replace adjacent 0s or 1s with 2s
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if input_grid[i][j] == 0 or input_grid[i][j] == 1:
                if (input_grid[i-1][j] == num and input_grid[i+1][j] == num) or \
                   (input_grid[i][j-1] == num and input_grid[i][j+1] == num):
                    output_grid[i][j] = 2

    return output_grid

# Test input grid
input_grid = [
    [5, 0, 0, 0, 0, 0, 5, 5],
    [5, 0, 5, 5, 5, 0, 5, 5],
    [5, 0, 0, 0, 0, 0, 5, 5],
    [5, 0, 5, 5, 5, 0, 5, 5],
    [5, 0, 0, 0, 0, 0, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))