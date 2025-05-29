def transform_grid(input_grid):
    # Determine the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Iterate over the grid to find blocks of different numbers
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check if the current cell is part of a block of different numbers
            if (input_grid[i][j] != input_grid[i-1][j] and
                input_grid[i][j] != input_grid[i+1][j] and
                input_grid[i][j] != input_grid[i][j-1] and
                input_grid[i][j] != input_grid[i][j+1]):
                # Replace the block with the surrounding pattern
                output_grid[i][j] = input_grid[i-1][j]
                output_grid[i][j+1] = input_grid[i-1][j+1]
                output_grid[i+1][j] = input_grid[i+1][j-1]
                output_grid[i+1][j+1] = input_grid[i+1][j]
    
    return output_grid

# Test input grid
input_grid = [
    [3, 0, 3, 0, 3, 0, 3, 0],
    [0, 3, 0, 3, 0, 3, 0, 3],
    [3, 0, 3, 0, 3, 0, 3, 0],
    [0, 3, 0, 3, 0, 3, 5, 5],
    [3, 0, 3, 0, 3, 0, 5, 5],
    [0, 3, 0, 3, 0, 3, 0, 3],
    [3, 0, 3, 0, 3, 0, 3, 0],
    [0, 3, 0, 3, 0, 3, 0, 3]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))