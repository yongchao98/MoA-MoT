input_grid = [0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Find the first and last non-zero elements
first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
last_non_zero_index = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1

# Calculate the number of leading zeros
leading_zeros = first_non_zero_index

# Calculate the number of positions to shift (leading zeros + 1)
shift_amount = leading_zeros + 1

# Create the output grid
output_grid = [0] * shift_amount + input_grid[first_non_zero_index:last_non_zero_index + 1] + [0] * (len(input_grid) - last_non_zero_index - shift_amount - 1)

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")