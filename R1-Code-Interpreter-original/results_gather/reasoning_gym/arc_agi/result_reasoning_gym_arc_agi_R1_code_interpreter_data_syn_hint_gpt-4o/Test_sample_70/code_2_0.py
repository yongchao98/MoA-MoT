# Define the input grid
input_grid = [
    [9, 9, 5],
    [7, 7, 4],
    [7, 7, 7]
]

# Initialize the output grid
output_grid = []

# Apply the transformation rule
# First row
output_grid.append(input_grid[0] + input_grid[0][::-1])

# Second row
output_grid.append(input_grid[1] + input_grid[1][::-1])

# Third and fourth rows
output_grid.append(input_grid[2] + input_grid[2])
output_grid.append(input_grid[2] + input_grid[2])

# Fifth row
output_grid.append(input_grid[1][::-1] + input_grid[1])

# Sixth row
output_grid.append(input_grid[0][::-1] + input_grid[0])

# Print the output grid
for row in output_grid:
    print(row)