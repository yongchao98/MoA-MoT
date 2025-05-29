def find_pattern(input_grid):
    # Convert input string to list
    grid = input_grid.split()
    grid = [int(x) for x in grid]
    
    # Find position of 2
    marker_pos = grid.index(2)
    
    # Create output grid
    output = grid.copy()
    
    # Check left side
    left_num = grid[marker_pos - 1] if marker_pos > 0 else 0
    right_num = grid[marker_pos + 1] if marker_pos < len(grid) - 1 else 0
    
    # If there's a sequence on the left
    if left_num != 0:
        output[marker_pos] = left_num
    # If there's a sequence on the right
    elif right_num != 0:
        output[marker_pos] = right_num
    
    # Convert back to string
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 2 0 3 3 3 3 3 3 3 3 3 3 3 0 0"
print(find_pattern(test_input))