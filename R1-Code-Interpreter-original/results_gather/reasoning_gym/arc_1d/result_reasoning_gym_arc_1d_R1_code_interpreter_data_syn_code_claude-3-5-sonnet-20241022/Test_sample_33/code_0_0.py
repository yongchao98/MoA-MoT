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
    
    # Find subsequent single numbers and expand them
    i = first_seq_len
    out_pos = first_seq_len
    while i < len(input_grid):
        if input_grid[i] != 0:
            # Found a number, expand it
            for j in range(first_seq_len):
                if out_pos + j < len(input_grid):
                    output[out_pos + j] = input_grid[i]
            out_pos += first_seq_len
        i += 1
    
    return output

# Test input
test_input = [2, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
result = find_pattern(test_input)
print(' '.join(map(str, result)))