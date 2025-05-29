# Define the input grid
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Find the first and last non-zero elements
first_non_zero_pos = next((i for i, x in enumerate(input_grid) if x != 0), None)
last_non_zero_pos = len(input_grid) - 1 - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None)

# Define the fixed shift
shift = 4

# Create the output grid
output_grid = [0] * len(input_grid)
for i in range(first_non_zero_pos, last_non_zero_pos + 1):
    new_pos = i + shift
    if new_pos < len(input_grid):
        output_grid[new_pos] = input_grid[i]

print(output_grid)