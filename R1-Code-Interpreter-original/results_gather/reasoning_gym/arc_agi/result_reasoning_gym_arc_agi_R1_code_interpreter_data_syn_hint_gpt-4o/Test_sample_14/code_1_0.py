def find_non_seven_cluster(input_grid):
    # Dimensions of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Iterate over the grid to find clusters of non-7 values
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check if the current position is part of a non-7 cluster
            if input_grid[i][j] != 7:
                # Extract the 3x3 grid centered at (i, j)
                output_grid = [
                    [input_grid[i-1][j-1], input_grid[i-1][j], input_grid[i-1][j+1]],
                    [input_grid[i][j-1], input_grid[i][j], input_grid[i][j+1]],
                    [input_grid[i+1][j-1], input_grid[i+1][j], input_grid[i+1][j+1]]
                ]
                return output_grid

# Test input grid
test_input = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 7, 7, 7, 7],
    [7, 7, 5, 7, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 5, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9],
    [7, 7, 5, 5, 5, 7, 7, 7, 7, 7, 7, 9, 9, 9],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 7, 7, 7, 9, 7],
    [7, 7, 7, 7, 7, 7, 7, 8, 7, 8, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 7, 7, 7, 7, 7],
    [7, 7, 7, 9, 9, 7, 7, 7, 7, 7, 7, 5, 7, 5],
    [7, 7, 9, 9, 9, 7, 7, 7, 7, 7, 7, 7, 5, 7],
    [7, 7, 7, 9, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

# Extract the output grid
output_grid = find_non_seven_cluster(test_input)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))