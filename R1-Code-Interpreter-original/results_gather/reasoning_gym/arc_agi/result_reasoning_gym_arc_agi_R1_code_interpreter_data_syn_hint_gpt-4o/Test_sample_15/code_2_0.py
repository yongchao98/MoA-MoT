# Define the input grid
input_grid = [
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 1, 2, 2, 2, 2],
    [2, 2, 2, 2, 0, 1, 0, 0, 0, 2],
    [2, 2, 2, 2, 1, 0, 1, 2, 2, 2],
    [2, 2, 2, 2, 0, 1, 0, 0, 0, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
]

# Initialize the output grid with the same values as the input grid
output_grid = [row[:] for row in input_grid]

# Grid dimensions
n = len(input_grid)
m = len(input_grid[0])

# Center lines
center_row = n // 2
center_col = m // 2

# Apply the mirroring rule
for i in range(n):
    for j in range(m):
        if input_grid[i][j] != 2:
            # Mirror across the vertical center line
            mirrored_col = m - 1 - j
            output_grid[i][mirrored_col] = input_grid[i][j]
            
            # Mirror across the horizontal center line
            mirrored_row = n - 1 - i
            output_grid[mirrored_row][j] = input_grid[i][j]
            
            # Mirror across both center lines
            output_grid[mirrored_row][mirrored_col] = input_grid[i][j]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))