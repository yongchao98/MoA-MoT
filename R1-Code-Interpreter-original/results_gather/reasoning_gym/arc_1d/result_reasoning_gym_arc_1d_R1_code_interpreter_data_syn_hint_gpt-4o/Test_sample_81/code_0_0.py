# Define the input grid
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Find the first non-zero element and its position
first_non_zero_pos = next((i for i, x in enumerate(input_grid) if x != 0), None)

# Calculate the shift
shift = len(input_grid) - first_non_zero_pos

# Create the output grid
output_grid = [0] * len(input_grid)
for i in range(first_non_zero_pos, len(input_grid)):
    if input_grid[i] != 0:
        new_pos = i + shift
        if new_pos < len(input_grid):
            output_grid[new_pos] = input_grid[i]

print(output_grid)