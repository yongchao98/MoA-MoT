def expand_grid(input_grid):
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Iterate over the grid to find blocks
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if input_grid[i][j] != 0:
                # Check if it's a central number
                if (input_grid[i-1][j] == input_grid[i+1][j] == input_grid[i][j-1] == input_grid[i][j+1] == input_grid[i][j]):
                    # Expand the block
                    num = input_grid[i][j]
                    outer_num = input_grid[i-1][j] if input_grid[i-1][j] != 0 else input_grid[i][j-1]
                    for di in [-1, 0, 1]:
                        for dj in [-1, 0, 1]:
                            if di != 0 or dj != 0:
                                output_grid[i+di][j+dj] = outer_num
                    output_grid[i][j] = num
    
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

output_grid = expand_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))