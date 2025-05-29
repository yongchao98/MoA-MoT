def expand_grid(input_grid):
    # Initialize an empty list for the output grid
    output_grid = []
    
    # Each element in the input grid is repeated to form a 4x4 grid
    for _ in range(4):
        output_grid.append([input_grid[0][0], input_grid[0][1], input_grid[0][0], input_grid[0][1]])
    
    return output_grid

# Define the test input grid
input_grid = [
    [4, 4],
    [4, 4]
]

# Get the output grid
output_grid = expand_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))