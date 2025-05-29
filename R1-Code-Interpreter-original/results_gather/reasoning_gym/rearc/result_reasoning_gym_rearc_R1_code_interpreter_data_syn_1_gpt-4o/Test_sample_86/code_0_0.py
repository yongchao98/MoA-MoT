# Define the input grid
input_grid = [
    [8, 4, 6, 8, 4],
    [8, 4, 6, 8, 4],
    [6, 6, 6, 6, 6],
    [8, 8, 6, 4, 8],
    [8, 4, 6, 4, 8]
]

# Identify the key row (all elements are the same)
key_row_index = 2  # The third row is the key row

# Construct the output grid
output_grid = [
    [1 if x == input_grid[key_row_index][0] else x for x in input_grid[key_row_index - 1]],
    [1 if x == input_grid[key_row_index][0] else x for x in input_grid[key_row_index + 1]]
]

# Print the output grid
print(output_grid)