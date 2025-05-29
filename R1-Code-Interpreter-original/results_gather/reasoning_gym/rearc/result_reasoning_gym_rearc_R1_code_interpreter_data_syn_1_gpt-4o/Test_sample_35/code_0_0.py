def transform_grid(input_grid):
    n = len(input_grid)
    m = len(input_grid[0])
    output_grid = [[0] * (2 * m) for _ in range(2 * n)]
    
    for i in range(n):
        for j in range(m):
            output_grid[i][j] = input_grid[i][j]
            output_grid[i][j + m] = input_grid[i][j]
            output_grid[i + n][j] = input_grid[i][j]
            output_grid[i + n][j + m] = input_grid[i][j]
    
    return output_grid

input_grid = [
    [0, 0],
    [0, 0]
]

output_grid = transform_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))