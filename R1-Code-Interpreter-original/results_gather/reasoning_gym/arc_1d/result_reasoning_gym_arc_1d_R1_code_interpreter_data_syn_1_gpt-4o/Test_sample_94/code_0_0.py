# Define the input grid
input_grid = [0, 0, 8, 3, 3, 1, 9, 2, 5, 7, 5, 2, 9, 2, 8, 3, 2, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Count the number of leading zeros
leading_zeros = 0
for num in input_grid:
    if num == 0:
        leading_zeros += 1
    else:
        break

# Count the number of non-zero elements after the initial zeros
non_zero_elements = input_grid[leading_zeros:]

# Calculate the number of leading zeros in the output
output_leading_zeros = leading_zeros + len(non_zero_elements)

# Create the output grid
output_grid = [0] * output_leading_zeros + non_zero_elements

# Print the output grid
print(output_grid)