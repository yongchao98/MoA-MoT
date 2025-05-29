from collections import Counter

def find_most_frequent_non_border(grid):
    # Flatten the grid and count occurrences of each number
    flat_grid = [num for row in grid for num in row]
    count = Counter(flat_grid)
    
    # Remove the border value (3) from the count
    if 3 in count:
        del count[3]
    
    # Find the most common non-border number
    most_common = count.most_common(1)
    return most_common[0][0] if most_common else None

# Test input grid
test_grid = [
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 7, 7, 7, 7, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 3, 3, 3, 3, 3, 3],
    [9, 9, 9, 9, 9, 9, 9, 9, 3, 2, 2, 2, 2, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3],
    [3, 6, 6, 6, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3],
    [3, 6, 6, 6, 5, 5, 5, 5, 3, 2, 2, 2, 2, 3, 3],
    [3, 6, 6, 6, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 5, 5, 5, 5, 3, 3, 4, 4, 4, 4, 3],
    [3, 3, 3, 3, 5, 5, 5, 5, 3, 3, 4, 3, 3, 4, 3],
    [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 3]
]

output = find_most_frequent_non_border(test_grid)
print(f"<<<{output}>>>")