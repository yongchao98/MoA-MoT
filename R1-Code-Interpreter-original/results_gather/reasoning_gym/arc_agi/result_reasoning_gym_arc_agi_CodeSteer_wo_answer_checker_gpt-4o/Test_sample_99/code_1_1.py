def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    # Function to check for cross pattern
    def is_cross_pattern(i, j, num):
        if i > 0 and i < rows - 1 and j > 0 and j < cols - 1:
            return (input_grid[i-1][j] == num and input_grid[i+1][j] == num and
                    input_grid[i][j-1] == num and input_grid[i][j+1] == num)
        return False

    # Function to check for diagonal pattern
    def is_diagonal_pattern(i, j, num):
        if i > 0 and i < rows - 1 and j > 0 and j < cols - 1:
            return (input_grid[i-1][j-1] == num and input_grid[i+1][j+1] == num and
                    input_grid[i-1][j+1] == num and input_grid[i+1][j-1] == num)
        return False

    # Apply transformation rules
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 5 and is_cross_pattern(i, j, 5):
                output_grid[i][j] = 9
            elif input_grid[i][j] == 4 and is_diagonal_pattern(i, j, 4):
                output_grid[i][j] = 2

    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 9, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [2, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 9, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))