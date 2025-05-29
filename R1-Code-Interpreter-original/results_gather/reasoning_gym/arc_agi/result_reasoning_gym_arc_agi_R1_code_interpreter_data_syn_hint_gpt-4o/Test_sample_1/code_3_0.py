def transform_grid(input_grid):
    n = len(input_grid)
    output_grid = [[0] * (2 * n) for _ in range(2 * n)]
    
    # Manually rearrange elements based on observed pattern
    # This is a placeholder for the actual pattern observed
    # The pattern needs to be identified and implemented here
    
    # Example pattern implementation (needs to be adjusted based on actual pattern)
    for i in range(n):
        for j in range(n):
            output_grid[i][j] = input_grid[j][i]
            output_grid[i][j + n] = input_grid[n - j - 1][i]
            output_grid[i + n][j] = input_grid[j][n - i - 1]
            output_grid[i + n][j + n] = input_grid[n - j - 1][n - i - 1]
    
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