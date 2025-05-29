# Test input
input_grid = [0, 1, 0, 0, 2, 1, 0, 0, 0, 1, 0, 1, 1, 2, 0]

# Count occurrences of 2, 1, and 0
count_2 = input_grid.count(2)
count_1 = input_grid.count(1)
count_0 = len(input_grid) - count_2 - count_1

# Construct the output grid
output_grid = [2] * count_2 + [1] * count_1 + [0] * count_0

# Print the output grid
print("<<<" + " ".join(map(str, output_grid)) + ">>>")