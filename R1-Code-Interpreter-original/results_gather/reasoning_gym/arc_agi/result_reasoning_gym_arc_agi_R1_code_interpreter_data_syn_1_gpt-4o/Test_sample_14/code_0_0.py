import numpy as np
from collections import Counter

def find_most_prominent_non_seven(grid):
    # Flatten the grid and count occurrences of each number
    flat_grid = grid.flatten()
    counter = Counter(flat_grid)
    # Remove 7 from the counter
    if 7 in counter:
        del counter[7]
    # Find the most common non-7 number
    most_common = counter.most_common(1)
    return most_common[0][0] if most_common else None

def extract_3x3_centered(grid, number):
    # Find all positions of the number in the grid
    positions = np.argwhere(grid == number)
    # Calculate the center position
    center_pos = positions[len(positions) // 2]
    # Extract the 3x3 grid centered around the center position
    r, c = center_pos
    return grid[r-1:r+2, c-1:c+2]

# Define the input grid
input_grid = np.array([
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 9, 9, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 7, 7, 7, 7],
    [7, 7, 5, 7, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 5, 7, 7, 7, 7, 7, 7, 7, 7, 9, 9],
    [7, 7, 5, 5, 5, 7, 7, 7, 7, 7, 7, 9, 9, 9],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 7, 7, 7, 9, 7],
    [7, 7, 7, 7, 7, 7, 7, 8, 7, 8, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 8, 7, 7, 7, 7, 7],
    [7, 7, 7, 9, 9, 7, 7, 7, 7, 7, 7, 5, 7, 5],
    [7, 7, 9, 9, 9, 7, 7, 7, 7, 7, 7, 7, 5, 7],
    [7, 7, 7, 9, 7, 7, 7, 7, 7, 7, 7, 5, 5, 5],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
])

# Find the most prominent non-7 number
prominent_number = find_most_prominent_non_seven(input_grid)

# Extract the 3x3 grid centered around the prominent number
output_grid = extract_3x3_centered(input_grid, prominent_number)

# Print the output grid
print(output_grid)