def transform_grid(input_grid):
    m, n = len(input_grid), len(input_grid[0])
    output_grid = [[0] * (2 * n) for _ in range(2 * m)]
    
    for i in range(m):
        for j in range(n):
            output_grid[2 * i][2 * j] = input_grid[i][j]
            output_grid[2 * i][2 * j + 1] = input_grid[i][j]
            output_grid[2 * i + 1][2 * j] = input_grid[i][j]
            output_grid[2 * i + 1][2 * j + 1] = input_grid[i][j]
    
    return output_grid

input_grid = [
    [9, 9, 3, 3],
    [9, 3, 9, 9],
    [9, 9, 9, 3],
    [3, 9, 3, 9]
]

output_grid = transform_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))