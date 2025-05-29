# Test input
input_grid = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Initialize the output grid with zeros
output_grid = [0] * len(input_grid)

# Apply the transformation rule
# Repeat the number 2 in the first few positions
output_grid[0:4] = [2, 2, 2, 2]

# Keep the number 1 in its position
output_grid[3:5] = [1, 1]

# Print the output grid
print(output_grid)