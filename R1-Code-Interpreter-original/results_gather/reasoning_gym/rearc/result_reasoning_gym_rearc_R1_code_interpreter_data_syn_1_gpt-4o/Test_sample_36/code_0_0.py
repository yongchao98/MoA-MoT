from collections import Counter

def find_least_frequent_number(grid):
    # Flatten the grid to a single list
    flat_grid = [num for row in grid for num in row]
    # Count the frequency of each number
    frequency = Counter(flat_grid)
    # Find the number with the least frequency that appears at least once
    least_frequent_number = min(frequency, key=lambda x: (frequency[x], x))
    return least_frequent_number

# Test input grid
test_input = [
    [0, 0, 0, 0, 0, 0, 0, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 0, 0, 0, 0, 0, 0, 0],
    [0, 7, 7, 7, 0, 0, 0, 0, 7],
    [0, 7, 6, 7, 0, 0, 0, 0, 0],
    [0, 7, 7, 7, 0, 7, 0, 0, 0],
    [6, 0, 0, 6, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Find and print the output for the test input
output = find_least_frequent_number(test_input)
print(output)