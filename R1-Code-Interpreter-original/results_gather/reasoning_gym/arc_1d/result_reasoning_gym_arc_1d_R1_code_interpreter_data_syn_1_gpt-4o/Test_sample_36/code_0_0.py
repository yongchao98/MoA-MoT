# Define the input grid
input_grid = [
    0, 0, 0, 9, 9,
    9, 9, 0, 0, 9,
    9, 9, 9, 9, 9,
    0, 0, 0, 9, 9,
    9, 9, 0, 0, 0
]

# Transform the third row to 1s
output_grid = input_grid[:10] + [1, 1, 1, 1, 1] + input_grid[15:]

# Print the output grid
print(output_grid)