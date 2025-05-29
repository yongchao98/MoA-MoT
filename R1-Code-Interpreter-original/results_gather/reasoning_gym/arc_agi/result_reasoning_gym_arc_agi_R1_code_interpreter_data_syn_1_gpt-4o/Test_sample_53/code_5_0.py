import numpy as np
from collections import Counter

# Test input grid
input_grid = np.array([
    [9, 4, 3, 7, 0, 2, 3, 3, 5, 7, 4, 1, 7, 8, 5, 7, 9, 8, 6],
    [8, 4, 8, 5, 9, 9, 4, 7, 4, 8, 0, 8, 7, 1, 9, 6, 6, 6, 9],
    [9, 3, 5, 1, 5, 7, 4, 9, 0, 8, 1, 8, 4, 2, 9, 5, 9, 5, 2],
    [0, 8, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 8, 6, 2, 7, 4, 4, 8],
    [9, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 4, 4, 1, 5, 9, 3],
    [7, 8, 3, 3, 3, 3, 3, 6, 3, 3, 3, 3, 3, 6, 3, 6, 9, 3, 8],
    [1, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 7, 1, 7, 6, 4, 6],
    [0, 6, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 7, 2, 1, 7],
    [1, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 5, 5, 3, 6, 7, 2],
    [5, 5, 3, 3, 6, 3, 3, 3, 3, 3, 3, 3, 5, 7, 5, 7, 8, 9, 3],
    [6, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 9, 8, 2, 0, 0, 3, 0],
    [1, 9, 8, 8, 5, 3, 3, 7, 7, 4, 4, 4, 7, 8, 5, 6, 8, 0, 5],
    [2, 0, 7, 8, 2, 1, 8, 1, 4, 8, 9, 3, 3, 5, 7, 1, 5, 3, 9],
    [1, 3, 6, 3, 6, 7, 6, 6, 0, 6, 4, 8, 8, 9, 6, 4, 1, 2, 3],
    [0, 2, 2, 6, 5, 8, 8, 6, 7, 2, 1, 8, 9, 1, 4, 3, 3, 1, 5],
    [4, 1, 9, 1, 4, 5, 1, 9, 8, 3, 4, 6, 3, 0, 7, 8, 9, 2, 6],
    [2, 9, 6, 7, 6, 2, 8, 5, 7, 4, 2, 3, 3, 9, 5, 8, 6, 5, 2]
])

# Flatten the grid and count occurrences of each number
flat_grid = input_grid.flatten()
counter = Counter(flat_grid)

# Find the two most common numbers
most_common = counter.most_common(2)
num1, num2 = most_common[0][0], most_common[1][0]

# Determine the size of the output grid
output_rows = 8
output_cols = 6

# Create the output grid
output_grid = np.full((output_rows, output_cols), num1)

# Fill the middle row with the second most common number
middle_row = output_rows // 2
output_grid[middle_row] = num2

# Alternate the pattern
for i in range(output_rows):
    if i % 2 == 0:
        output_grid[i] = num1
    else:
        output_grid[i] = num2

# Print the output grid
print(output_grid)