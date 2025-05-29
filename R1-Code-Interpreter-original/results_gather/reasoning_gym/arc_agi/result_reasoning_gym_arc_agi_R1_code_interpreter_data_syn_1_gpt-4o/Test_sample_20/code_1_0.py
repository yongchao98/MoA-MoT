def expand_grid(input_grid):
    # Determine the size of the input grid
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # Initialize the output grid with the appropriate size
    output_rows = input_rows * 4
    output_cols = input_cols * 4
    output_grid = [[2 for _ in range(output_cols)] for _ in range(output_rows)]
    
    # Fill the output grid based on the input grid
    for i in range(input_rows):
        for j in range(input_cols):
            # The center of the 4x4 block in the output grid
            for di in range(4):
                for dj in range(4):
                    output_grid[i*4 + di][j*4 + dj] = 2
            output_grid[i*4 + 1][j*4 + 1] = input_grid[i][j]
            output_grid[i*4 + 1][j*4 + 2] = input_grid[i][j]
            output_grid[i*4 + 2][j*4 + 1] = input_grid[i][j]
            output_grid[i*4 + 2][j*4 + 2] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [3, 2, 2, 3, 2, 5, 5, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [4, 2, 2, 3, 2, 5, 5, 5]
]

# Generate the output grid
output_grid = expand_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))