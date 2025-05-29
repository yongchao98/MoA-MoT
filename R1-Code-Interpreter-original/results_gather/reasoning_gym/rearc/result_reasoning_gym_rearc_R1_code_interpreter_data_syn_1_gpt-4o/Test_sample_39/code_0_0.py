def find_output_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Iterate over the grid to find the distinct number
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check if the current number is distinct
            if input_grid[i][j] != input_grid[i-1][j] and input_grid[i][j] != input_grid[i+1][j] and \
               input_grid[i][j] != input_grid[i][j-1] and input_grid[i][j] != input_grid[i][j+1]:
                # Extract the 3x3 grid centered around the distinct number
                output_grid = [
                    [input_grid[i-1][j-1], input_grid[i-1][j], input_grid[i-1][j+1]],
                    [input_grid[i][j-1], input_grid[i][j], input_grid[i][j+1]],
                    [input_grid[i+1][j-1], input_grid[i+1][j], input_grid[i+1][j+1]]
                ]
                return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4],
    [4, 4, 4, 7, 7, 7],
    [4, 4, 4, 7, 2, 7],
    [4, 4, 4, 7, 7, 7]
]

output_grid = find_output_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))