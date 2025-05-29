from collections import Counter

def find_least_frequent_count(grid):
    # Flatten the grid and count occurrences of each number
    flat_grid = [num for row in grid for num in row]
    counts = Counter(flat_grid)
    # Find the minimum count
    min_count = min(counts.values())
    return min_count

# Test input grid
test_input = [
    [4, 3, 3],
    [4, 3, 3],
    [3, 4, 4]
]

# Find the output for the test input
output = find_least_frequent_count(test_input)
print(output)