def find_most_frequent_number(grid):
    from collections import Counter
    
    # Flatten the grid and filter out zeros
    numbers = [num for row in grid for num in row if num != 0]
    
    # Count the frequency of each number
    frequency = Counter(numbers)
    
    # Find the number with the highest frequency
    most_frequent_number = frequency.most_common(1)[0][0]
    
    return most_frequent_number

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 0, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7]
]

# Find and print the output
output = find_most_frequent_number(input_grid)
print(output)