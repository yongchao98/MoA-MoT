# Test input grid
input_grid = [2, 2, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 0, 0, 0, 0]

# Count occurrences of 2s, 0s, and 1s
count_2 = input_grid.count(2)
count_0 = input_grid.count(0)
count_1 = input_grid.count(1)

# Construct the output grid
output_grid = [2] * count_2 + [0] * count_0 + [1] * count_1

# Print the output grid
print(output_grid)