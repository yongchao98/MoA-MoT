def expand_grid(input_grid):
    # Determine the size of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Initialize an empty list for the output grid
    # The output grid will be larger to accommodate the 3x3 blocks
    output_grid = [[0] * cols for _ in range(rows)]
    
    # Iterate over each cell in the input grid
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 0:
                # Get the number to expand
                num = input_grid[i][j]
                
                # Place the 3x3 block in the output grid
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        ni, nj = i + di, j + dj
                        if 0 <= ni < rows and 0 <= nj < cols:
                            output_grid[ni][nj] = num
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 5, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Get the output grid
output_grid = expand_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))