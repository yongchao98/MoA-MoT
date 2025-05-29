def shift_grid(input_grid):
    # Create an empty grid for the output
    output_grid = [[4 for _ in range(10)] for _ in range(10)]
    
    # Define the shift pattern observed from examples
    shift_down = 2
    shift_right = 2
    
    # Iterate over the input grid to find clusters and shift them
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 4:
                new_i = (i + shift_down) % 10
                new_j = (j + shift_right) % 10
                output_grid[new_i][new_j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 8, 8, 4, 4, 4, 8, 8, 4],
    [4, 8, 8, 8, 8, 4, 8, 8, 8, 8],
    [4, 4, 8, 4, 4, 4, 4, 8, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 2, 2, 4],
    [4, 4, 8, 8, 4, 4, 1, 1, 2, 1],
    [4, 8, 8, 8, 8, 4, 4, 1, 4, 4],
    [4, 4, 8, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
]

# Get the output grid
output_grid = shift_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))