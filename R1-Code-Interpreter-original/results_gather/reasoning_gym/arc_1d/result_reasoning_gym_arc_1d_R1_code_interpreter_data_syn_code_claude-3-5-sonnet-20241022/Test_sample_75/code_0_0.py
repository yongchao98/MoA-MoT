def find_pattern(input_grid):
    # Convert input to list
    grid = [int(x) for x in input_grid.split()]
    
    # Find first non-zero sequence from left
    start_idx = 0
    while start_idx < len(grid) and grid[start_idx] == 0:
        start_idx += 1
    
    # Find length of non-zero sequence
    seq_len = 0
    for i in range(start_idx, len(grid)):
        if grid[i] != 0:
            seq_len += 1
        else:
            break
            
    # If all zeros or no clear sequence found
    if seq_len == 0:
        return grid
    
    # Create output by circular shift
    result = grid[seq_len:] + grid[:seq_len]
    
    # Convert to string
    return ' '.join(map(str, result))

# Test input
test_input = "5 5 9 7 4 1 5 6 9 2 7 8 1 7 1 6 8 9 9 8 8 5 4 7 4 9 6 7 7"
print(find_pattern(test_input))