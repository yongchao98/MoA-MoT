def transform_grid(input_grid):
    n = len(input_grid)
    output_grid = [[0] * (2 * n) for _ in range(2 * n)]
    
    for i in range(n):
        for j in range(n):
            # Create a 2x2 block for each element
            output_grid[2 * i][2 * j] = input_grid[i][j]
            output_grid[2 * i][2 * j + 1] = input_grid[i][j]
            output_grid[2 * i + 1][2 * j] = input_grid[i][j]
            output_grid[2 * i + 1][2 * j + 1] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [4, 6, 4, 2],
    [2, 5, 2, 4],
    [6, 5, 6, 6],
    [6, 5, 6, 6]
]

# Generate the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))