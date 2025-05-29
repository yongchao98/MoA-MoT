def transform_grid(input_grid):
    output_grid = [[0] * 10 for _ in range(10)]
    
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 2:
                output_grid[i][j] = 2
            elif input_grid[i][j] == 8:
                # Apply the transformation rule for 8
                new_i = i + 2
                new_j = j + 2
                if new_i < 10 and new_j < 10:
                    output_grid[new_i][new_j] = 8
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 8, 0],
    [0, 0, 0, 8, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 2, 0, 0, 0, 0, 0, 0],
    [8, 0, 2, 2, 0, 0, 0, 0, 0, 8],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))