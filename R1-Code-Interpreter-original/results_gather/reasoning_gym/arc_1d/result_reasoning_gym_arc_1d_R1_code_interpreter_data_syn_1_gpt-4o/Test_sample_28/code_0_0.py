# Test input grid
input_grid = [1, 0, 0, 1, 2, 1, 0, 0, 0, 0, 1, 1, 2, 0, 2, 0, 0, 2, 0]

# Count the number of 1s and 2s
count_1 = input_grid.count(1)
count_2 = input_grid.count(2)

# Calculate the number of 0s
count_0 = len(input_grid) - count_1 - count_2

# Construct the output grid
output_grid = [1] * count_1 + [0] * count_0 + [2] * count_2

# Print the output grid
print(output_grid)