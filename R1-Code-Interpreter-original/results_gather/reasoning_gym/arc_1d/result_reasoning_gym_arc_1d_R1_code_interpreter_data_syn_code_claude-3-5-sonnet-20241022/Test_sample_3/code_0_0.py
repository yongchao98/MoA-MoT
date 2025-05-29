def find_pattern(input_grid):
    # Convert to list for manipulation
    grid = list(map(int, input_grid.split()))
    
    # Find position of 1
    one_pos = grid.index(1)
    
    # Count non-zero numbers after 1
    non_zeros = [x for x in grid[one_pos+1:] if x != 0]
    
    # Create output grid
    output = [0] * len(grid)
    output[one_pos] = 1  # Keep 1 in same position
    
    # Place non-zeros right after 1
    for i, val in enumerate(non_zeros):
        if one_pos + 1 + i < len(grid):
            output[one_pos + 1 + i] = val
    
    # Convert back to string
    return ' '.join(map(str, output))

# Test input
test_input = "1 0 0 0 9 0 0 0 0 0 0 0 0 0"
print(find_pattern(test_input))