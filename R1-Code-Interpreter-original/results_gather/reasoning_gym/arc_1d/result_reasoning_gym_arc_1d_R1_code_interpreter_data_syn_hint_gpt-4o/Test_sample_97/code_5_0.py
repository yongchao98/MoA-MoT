# Test input
input_grid = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Initialize the output grid with zeros
output_grid = [0] * len(input_grid)

# Find the first non-zero number and its position
first_non_zero = None
first_non_zero_pos = None

for i, num in enumerate(input_grid):
    if num != 0:
        first_non_zero = num
        first_non_zero_pos = i
        break

# Repeat the first non-zero number up to its position
if first_non_zero is not None:
    output_grid[:first_non_zero_pos] = [first_non_zero] * first_non_zero_pos

# Copy the rest of the input grid to the output grid
output_grid[first_non_zero_pos:] = input_grid[first_non_zero_pos:]

# Print the output grid
print(output_grid)