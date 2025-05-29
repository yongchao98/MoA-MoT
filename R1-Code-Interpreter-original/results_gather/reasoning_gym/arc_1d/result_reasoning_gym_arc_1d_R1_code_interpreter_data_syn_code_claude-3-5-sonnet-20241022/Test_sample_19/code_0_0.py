def find_pattern(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Find first non-zero sequence
    start = 0
    while start < len(grid) and grid[start] == 0:
        start += 1
    
    if start == len(grid):
        return grid  # All zeros
    
    # Find end of non-zero sequence
    end = start
    while end < len(grid) and (grid[end] != 0 or (end > 0 and grid[end-1] != 0)):
        end += 1
    
    # Extract sequence
    sequence = grid[start:end]
    
    # Create output grid
    output = sequence + [0] * (len(grid) - len(sequence))
    
    # Convert to string format
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 0 0 1 8 9 2 8 8 6 3 9 8 2 6 4 7 7"
print(find_pattern(test_input))