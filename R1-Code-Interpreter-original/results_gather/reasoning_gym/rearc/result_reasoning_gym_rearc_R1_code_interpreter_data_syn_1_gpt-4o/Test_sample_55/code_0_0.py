def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Identify the block of different numbers
    block_value = None
    for row in input_grid:
        for value in row:
            if input_grid[0][0] != value:
                block_value = value
                break
        if block_value is not None:
            break
    
    # Change the adjacent region to 4
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] == block_value:
                # Check adjacent cells and change them to 4 if they are not the block value
                if i > 0 and input_grid[i-1][j] != block_value:
                    output_grid[i-1][j] = 4
                if i < len(input_grid) - 1 and input_grid[i+1][j] != block_value:
                    output_grid[i+1][j] = 4
                if j > 0 and input_grid[i][j-1] != block_value:
                    output_grid[i][j-1] = 4
                if j < len(input_grid[i]) - 1 and input_grid[i][j+1] != block_value:
                    output_grid[i][j+1] = 4
    
    return output_grid

# Test input grid
input_grid = [
    [2, 2, 2, 2, 2, 2, 3, 2],
    [2, 2, 2, 2, 2, 2, 3, 2],
    [2, 2, 2, 2, 3, 3, 3, 3],
    [2, 2, 2, 2, 3, 2, 2, 3],
    [2, 2, 2, 2, 3, 2, 2, 3],
    [3, 3, 3, 3, 3, 2, 2, 3],
    [2, 2, 2, 2, 3, 2, 2, 3],
    [2, 2, 2, 2, 3, 3, 3, 3]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))