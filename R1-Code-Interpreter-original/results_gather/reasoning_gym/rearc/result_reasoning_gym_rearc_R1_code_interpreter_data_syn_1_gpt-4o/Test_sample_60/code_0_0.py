# Define the input grid
input_grid = [
    [1, 1, 8, 8, 8, 1],
    [1, 1, 8, 7, 8, 1],
    [1, 1, 8, 8, 8, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1]
]

# Extract the center number
center_number = input_grid[1][3]

# Determine the most frequent number around the center
from collections import Counter

# Get the surrounding numbers
surrounding_numbers = [
    input_grid[0][2], input_grid[0][3], input_grid[0][4],
    input_grid[1][2], input_grid[1][4],
    input_grid[2][2], input_grid[2][3], input_grid[2][4]
]

# Find the most common surrounding number
most_common_number = Counter(surrounding_numbers).most_common(1)[0][0]

# Create the output grid
output_grid = [
    [most_common_number, most_common_number, most_common_number],
    [most_common_number, center_number, most_common_number],
    [most_common_number, most_common_number, most_common_number]
]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))