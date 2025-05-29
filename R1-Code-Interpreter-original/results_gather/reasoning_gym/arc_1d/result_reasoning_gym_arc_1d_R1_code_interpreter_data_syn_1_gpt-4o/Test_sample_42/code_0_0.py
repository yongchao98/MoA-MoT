# Test input grid
input_grid = [9, 0, 9, 9, 9, 0, 9, 9, 0, 0, 0, 0, 0]

# Applying the derived rule
output_grid = [
    input_grid[0],  # Output[0]
    0,              # Output[1]
    input_grid[2],  # Output[2]
    input_grid[0],  # Output[3]
    0,              # Output[4]
    input_grid[2],  # Output[5]
    input_grid[0],  # Output[6]
    input_grid[2],  # Output[7]
    0, 0, 0, 0, 0   # Output[8] to Output[12]
]

print(output_grid)