def transform_grid(input_grid):
    # Create an empty grid for the output
    output_grid = [[0] * 10 for _ in range(10)]
    
    # Copy '2's to the output grid
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 2:
                output_grid[i][j] = 2
    
    # Move '8's based on observed pattern
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] == 8:
                # Determine new position for '8'
                # Example pattern: move diagonally down-right by 2 positions if possible
                new_i, new_j = i + 2, j + 2
                if new_i < 10 and new_j < 10 and output_grid[new_i][new_j] == 0:
                    output_grid[new_i][new_j] = 8
                else:
                    # If not possible, try moving down by 2
                    new_i, new_j = i + 2, j
                    if new_i < 10 and output_grid[new_i][new_j] == 0:
                        output_grid[new_i][new_j] = 8
                    else:
                        # If not possible, try moving right by 2
                        new_i, new_j = i, j + 2
                        if new_j < 10 and output_grid[new_i][new_j] == 0:
                            output_grid[new_i][new_j] = 8
    
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