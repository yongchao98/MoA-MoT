def analyze_pattern(input_grid, output_grid):
    # Count leading zeros
    input_leading = 0
    output_leading = 0
    for i in range(len(input_grid)):
        if input_grid[i] == 0:
            input_leading += 1
        else:
            break
    for i in range(len(output_grid)):
        if output_grid[i] == 0:
            output_leading += 1
        else:
            break
            
    # Count sequence length
    seq_length = 0
    for i in range(input_leading, len(input_grid)):
        if input_grid[i] != 0:
            seq_length += 1
        else:
            break
            
    # Apply pattern to test input
    test_input = [0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    # Count leading zeros and sequence in test input
    test_leading = 0
    test_seq = 0
    for i in range(len(test_input)):
        if test_input[i] == 0:
            test_leading += 1
        else:
            break
    for i in range(test_leading, len(test_input)):
        if test_input[i] != 0:
            test_seq += 1
        else:
            break
    
    # Generate output based on pattern
    result = []
    # Add leading zeros (double the original leading zeros)
    result.extend([0] * (test_leading * 2))
    # Add the sequence
    result.extend([8] * test_seq)
    # Fill remaining with zeros
    while len(result) < 23:  # Leave space for final zero
        result.append(0)
    # Add final zero
    result.append(0)
    
    print(' '.join(map(str, result)))

# Test examples
example1_in = [0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0]
example1_out = [0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0]

analyze_pattern(example1_in, example1_out)