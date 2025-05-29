# Test input grid
input_grid = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0]

# Initialize the output grid with zeros
output_grid = [0] * len(input_grid)

# Set the first and last non-zero elements
output_grid[0] = input_grid[0]
output_grid[17] = input_grid[17]

# Print the output grid
print(output_grid)