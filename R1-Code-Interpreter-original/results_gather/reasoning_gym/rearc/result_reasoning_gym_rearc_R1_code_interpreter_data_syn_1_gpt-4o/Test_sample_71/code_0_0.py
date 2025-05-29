from collections import Counter

# Test input grid
input_grid = [
    [5, 5, 5, 2, 5, 5, 2, 5, 5, 5],
    [5, 5, 5, 2, 2, 2, 2, 5, 5, 5],
    [5, 5, 5, 0, 5, 5, 0, 5, 5, 5],
    [2, 2, 0, 5, 0, 0, 5, 0, 2, 2],
    [5, 2, 5, 0, 5, 5, 0, 5, 2, 5],
    [5, 2, 5, 0, 5, 5, 0, 5, 2, 5],
    [2, 2, 0, 5, 0, 0, 5, 0, 2, 2],
    [5, 5, 5, 3, 5, 5, 0, 5, 5, 5],
    [5, 5, 5, 2, 2, 2, 2, 5, 5, 5],
    [5, 5, 5, 2, 5, 5, 2, 5, 5, 5]
]

# Flatten the grid to a single list
flat_list = [item for sublist in input_grid for item in sublist]

# Count the frequency of each number
counter = Counter(flat_list)

# Find the most common number
most_common_number = counter.most_common(1)[0][0]

print(most_common_number)