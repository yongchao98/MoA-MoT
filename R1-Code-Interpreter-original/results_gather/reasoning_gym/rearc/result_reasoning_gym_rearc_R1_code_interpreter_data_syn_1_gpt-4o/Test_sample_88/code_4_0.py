def transform_grid(input_grid):
    output_grid = [row[:] for row in input_grid]  # Copy the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 7:  # Only transform numbers that are not 7
                if i == rows - 1 or j == cols - 1:  # Bottom row or rightmost column
                    output_grid[i][j] = 3
                else:
                    output_grid[i][j] = 2
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 8, 8, 7],
    [7, 7, 7, 8, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 8, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 8, 8, 7, 7, 7, 7, 7],
    [8, 7, 7, 7, 7, 7, 7, 7, 7],
    [8, 7, 7, 7, 7, 7, 7, 7, 7],
    [8, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))