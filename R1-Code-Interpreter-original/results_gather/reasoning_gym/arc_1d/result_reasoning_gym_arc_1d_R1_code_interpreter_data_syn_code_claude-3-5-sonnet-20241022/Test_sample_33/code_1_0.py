def find_pattern(input_grid):
    # Find the length of the first sequence
    first_seq_len = 0
    num = input_grid[0]
    for i in range(len(input_grid)):
        if input_grid[i] != num:
            break
        first_seq_len += 1
    
    # Create output grid
    output = [0] * len(input_grid)
    
    # Copy the first sequence
    for i in range(first_seq_len):
        output[i] = input_grid[i]
    
    # Process the rest
    current_pos = first_seq_len
    while current_pos < len(input_grid):
        if input_grid[current_pos] != 0:
            # Found a number, add zeros until next position
            next_pos = current_pos
            # Create a group of length first_seq_len
            for i in range(first_seq_len):
                if next_pos < len(input_grid):
                    output[next_pos] = input_grid[current_pos]
                    next_pos += 1
        current_pos += 1
    
    return output

# Test input
test_input = [2, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
result = find_pattern(test_input)
print(' '.join(map(str, result)))