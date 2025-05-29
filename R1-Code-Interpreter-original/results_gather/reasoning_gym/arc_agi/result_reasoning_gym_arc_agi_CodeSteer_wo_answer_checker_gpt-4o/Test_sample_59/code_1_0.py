def transform_grid(input_grid):
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Define the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the central block
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if input_grid[i][j] != 0:
                # Check if it's the start of a central block
                if (input_grid[i-1][j] == 0 and input_grid[i][j-1] == 0 and
                    input_grid[i+1][j] != 0 and input_grid[i][j+1] != 0):
                    
                    # Determine the size of the central block
                    central_value = input_grid[i][j]
                    size = 0
                    while i + size < rows and input_grid[i + size][j] == central_value:
                        size += 1
                    
                    # Expand the central block
                    for x in range(i - 1, i + size + 1):
                        for y in range(j - 1, j + size + 1):
                            output_grid[x][y] = central_value
                    
                    # Find and expand the smaller block
                    smaller_value = input_grid[i + 1][j + 1]
                    for x in range(i, i + size):
                        for y in range(j, j + size):
                            if input_grid[x][y] == smaller_value:
                                output_grid[x][y] = smaller_value
                                output_grid[x-1][y] = smaller_value
                                output_grid[x+1][y] = smaller_value
                                output_grid[x][y-1] = smaller_value
                                output_grid[x][y+1] = smaller_value
                    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 4, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 4, 4, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 4, 4, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))