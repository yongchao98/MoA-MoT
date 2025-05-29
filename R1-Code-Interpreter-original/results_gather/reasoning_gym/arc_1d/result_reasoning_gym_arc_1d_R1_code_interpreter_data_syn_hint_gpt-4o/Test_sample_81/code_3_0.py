# Define the input grid
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Find the first and last non-zero elements
first_non_zero_pos = next((i for i, x in enumerate(input_grid) if x != 0), None)
last_non_zero_pos = len(input_grid) - 1 - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)

# Calculate the length of the non-zero sequence
non_zero_length = last_non_zero_pos - first_non_zero_pos + 1

# Calculate the shift based on the position of the first non-zero element
# Hypothesis: The shift is calculated as (length of grid - position of first non-zero element - non-zero length)
shift = len(input_grid) - first_non_zero_pos - non_zero_length

# Create the output grid
output_grid = [0] * len(input_grid)
for i in range(first_non_zero_pos, last_non_zero_pos + 1):
    new_pos = i + shift
    if new_pos < len(input_grid):
        output_grid[new_pos] = input_grid[i]

print(output_grid)