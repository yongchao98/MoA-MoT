def find_pattern(input_grid):
    # Find the length of the first sequence
    first_seq_len = 0
    num = input_grid[0]
    for i in range(len(input_grid)):
        if input_grid[i] != num:
            break
        first_seq_len += 1
    
    output = [0] * len(input_grid)
    
    # Copy first sequence
    for i in range(first_seq_len):
        output[i] = input_grid[i]
    
    # Find positions of single numbers in input
    positions = []
    for i in range(first_seq_len, len(input_grid)):
        if input_grid[i] != 0:
            positions.append(i)
    
    # Process each position
    current_pos = first_seq_len
    for pos in positions:
        # Add one zero after previous group
        current_pos += 1
        # Add the new group
        for i in range(first_seq_len):
            if current_pos < len(input_grid):
                output[current_pos] = input_grid[pos]
                current_pos += 1
    
    return output

# Test input
test_input = [2, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
result = find_pattern(test_input)
print(' '.join(map(str, result)))