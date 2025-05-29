# Define the input grid
input_grid = [
    [0, 0, 3, 2, 3, 3, 3, 3],
    [0, 0, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 2, 2, 3],
    [3, 2, 2, 3, 3, 2, 2, 3],
    [3, 2, 2, 3, 0, 3, 3, 3]
]

# Extract the 3x3 section from the top-left corner
output_grid = [row[:3] for row in input_grid[:3]]

# Set the bottom-right corner to 9
output_grid[2][2] = 9

# Print the output grid
print(output_grid)