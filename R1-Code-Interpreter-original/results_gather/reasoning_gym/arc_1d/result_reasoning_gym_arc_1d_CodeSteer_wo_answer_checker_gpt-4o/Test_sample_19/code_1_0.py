input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 8, 9, 2, 8, 8, 6, 3, 9, 8, 2, 6, 4, 7, 7]

# Find the first non-zero sequence
start_index = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
sequence = input_grid[start_index:]

# Calculate the number of leading zeros
leading_zeros = start_index

# Create the output grid by shifting the sequence
output_grid = sequence + [0] * leading_zeros

# Replace the last zero with the first number of the sequence if there are zeros left
if leading_zeros > 0:
    output_grid[-1] = sequence[0]

print("<<<" + " ".join(map(str, output_grid)) + ">>>")