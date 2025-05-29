def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[7] * cols for _ in range(rows)]
    
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if input_grid[i][j] != 7:
                # Check if the current cell is surrounded by 7s
                if (input_grid[i-1][j] == 7 and input_grid[i+1][j] == 7 and
                    input_grid[i][j-1] == 7 and input_grid[i][j+1] == 7):
                    output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [6, 7, 7, 7, 6, 7, 7, 7, 7, 6],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 6, 7, 6, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [6, 7, 6, 6, 6, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 6, 7, 7, 7, 7, 6, 6, 7, 7],
    [7, 7, 7, 7, 7, 7, 6, 6, 6, 7],
    [7, 7, 6, 7, 7, 6, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))