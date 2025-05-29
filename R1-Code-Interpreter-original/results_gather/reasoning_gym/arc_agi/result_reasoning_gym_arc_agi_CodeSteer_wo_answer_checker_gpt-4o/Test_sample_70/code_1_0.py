# Define the input grid
input_grid = [
    [9, 9, 5],
    [7, 7, 4],
    [7, 7, 7]
]

# Create the output grid based on the identified pattern
output_grid = []

# Append the first row and its reverse
output_grid.append(input_grid[0] + input_grid[0][::-1])

# Append the second row and its reverse
output_grid.append(input_grid[1] + input_grid[1][::-1])

# Append the third row and its reverse
output_grid.append(input_grid[2] + input_grid[2][::-1])

# Repeat the third row
output_grid.append(input_grid[2] + input_grid[2][::-1])

# Repeat the second row
output_grid.append(input_grid[1] + input_grid[1][::-1])

# Repeat the first row
output_grid.append(input_grid[0] + input_grid[0][::-1])

# Print the output grid
for row in output_grid:
    print(row)