def find_pattern(input_grid):
    # Convert input string to list
    grid = [int(x) for x in input_grid.split()]
    output = grid.copy()
    
    # Find position of 2
    marker_pos = grid.index(2)
    
    # Look left and right for non-zero sequences
    left_sequence = None
    right_sequence = None
    
    # Check left side
    for i in range(marker_pos - 1, -1, -1):
        if grid[i] != 0:
            left_sequence = grid[i]
            break
            
    # Check right side
    for i in range(marker_pos + 1, len(grid)):
        if grid[i] != 0:
            right_sequence = grid[i]
            break
    
    # Apply the rule
    if right_sequence is not None:
        # Extend right sequence one position left
        output[marker_pos] = right_sequence
    elif left_sequence is not None:
        # Extend left sequence one position right
        output[marker_pos] = left_sequence
    
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 2 0 3 3 3 3 3 3 3 3 3 3 3 0 0"
print(find_pattern(test_input))