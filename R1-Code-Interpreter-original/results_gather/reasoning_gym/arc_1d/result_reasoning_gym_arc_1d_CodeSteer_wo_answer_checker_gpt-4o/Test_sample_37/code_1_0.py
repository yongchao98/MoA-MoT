# Test input
input_grid = [0, 0, 4, 3, 5, 4, 9, 0, 0, 0, 0, 0, 0]

# Find the first non-zero index
first_non_zero_index = next(i for i, x in enumerate(input_grid) if x != 0)

# Create the output grid
output_grid = input_grid[:]
output_grid[first_non_zero_index] = 0
output_grid[first_non_zero_index + 1:first_non_zero_index + 6] = input_grid[first_non_zero_index:first_non_zero_index + 5]

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")