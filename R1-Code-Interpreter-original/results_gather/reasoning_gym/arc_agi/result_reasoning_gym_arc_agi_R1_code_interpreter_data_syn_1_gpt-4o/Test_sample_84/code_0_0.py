# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 6, 6, 6, 9, 9, 9],
    [9, 9, 9, 6, 9, 6, 6, 6, 9],
    [9, 6, 6, 6, 9, 9, 9, 6, 6],
    [6, 6, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 6, 6, 6, 9, 9],
    [6, 6, 9, 9, 6, 9, 6, 6, 6],
    [9, 6, 6, 6, 6, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9]
]

# Determine the number of '9's in the first row
first_row = input_grid[0]
num_nines = first_row.count(9)

# Calculate the length of the output grid
output_length = num_nines // 2

# Generate the output grid
output_grid = [9] * output_length

# Print the output grid
print(output_grid)