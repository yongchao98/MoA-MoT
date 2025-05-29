def transform_grid(input_grid):
    # Define the size of the input and output grids
    n = 4
    m = 8
    
    # Initialize the output grid
    output_grid = [[0] * m for _ in range(m)]
    
    # Fill the output grid based on the observed pattern
    for i in range(n):
        for j in range(n):
            # Top-left quadrant
            output_grid[i][j] = input_grid[i][j]
            # Top-right quadrant
            output_grid[i][j + n] = input_grid[j][n - i - 1]
            # Bottom-left quadrant
            output_grid[i + n][j] = input_grid[n - j - 1][i]
            # Bottom-right quadrant
            output_grid[i + n][j + n] = input_grid[n - i - 1][n - j - 1]
    
    return output_grid

# Test the function with the provided test input
test_input = [
    [4, 6, 4, 2],
    [2, 5, 2, 4],
    [6, 5, 6, 6],
    [6, 5, 6, 6]
]

# Transform the test input grid
output_grid = transform_grid(test_input)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))