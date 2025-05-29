# Test input grid
input_grid = [0, 1, 0, 0, 2, 1, 0, 0, 0, 1, 0, 1, 1, 2, 0]

# Count the number of 2's and 1's
count_2 = input_grid.count(2)
count_1 = input_grid.count(1)

# Construct the output grid
output_grid = [2] * count_2 + [1] * count_1 + [0] * (len(input_grid) - count_2 - count_1)

# Print the output grid
print(output_grid)