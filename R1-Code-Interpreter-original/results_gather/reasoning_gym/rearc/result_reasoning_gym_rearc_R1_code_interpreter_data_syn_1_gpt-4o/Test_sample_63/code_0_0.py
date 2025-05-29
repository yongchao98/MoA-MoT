def find_unique_largest(grid):
    from collections import Counter
    
    # Flatten the grid and count the frequency of each number
    flat_grid = [num for row in grid for num in row]
    frequency = Counter(flat_grid)
    
    # Find the largest number that appears only once
    unique_largest = max((num for num, count in frequency.items() if count == 1), default=None)
    
    return unique_largest

# Test input grid
test_grid = [
    [6, 6, 6, 6],
    [6, 2, 6, 2],
    [6, 2, 6, 2],
    [6, 9, 6, 6]
]

# Find the output for the test input grid
output = find_unique_largest(test_grid)
print(output)