def transform_grid(input_grid):
    # Create an empty grid for the output
    output_grid = [[0] * 10 for _ in range(10)]
    
    # Copy '2's to the output grid
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 2:
                output_grid[i][j] = 2
    
    # Move '8's based on a more complex pattern
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 8:
                # Attempt to move '8' to a new position
                # Prioritize moving to a position that is two rows down or two columns to the right
                new_positions = [
                    (i + 3, j - 1),  # Move down and left
                    (i + 2, j + 1),  # Move down and right
                    (i + 1, j + 2),  # Move right and down
                    (i + 2, j - 2)   # Move down and left
                ]
                
                for new_i, new_j in new_positions:
                    if 0 <= new_i < 10 and 0 <= new_j < 10 and output_grid[new_i][new_j] == 0:
                        output_grid[new_i][new_j] = 8
                        break
    
    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 8, 0],
    [0, 0, 0, 8, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 2, 0, 0, 0, 0, 0, 0],
    [8, 0, 2, 2, 0, 0, 0, 0, 0, 8],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))