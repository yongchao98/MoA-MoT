# Test input grid
input_grid = [0, 0, 0, 2, 0, 2, 1, 0, 0, 0, 1, 1, 1, 0, 2, 1, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0]

# Count the occurrences of each number
count_zeros = input_grid.count(0)
count_ones = input_grid.count(1)
count_twos = input_grid.count(2)

# Construct the output grid based on the rule
output_grid = [0] * count_zeros + [1] * count_ones + [2] * count_twos

# Print the output grid
print(output_grid)