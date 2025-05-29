def transform_grid(input_grid):
    # Find position of 2
    two_pos = input_grid.index(2)
    
    # Find the sequence of non-zero numbers (excluding 2)
    sequence = []
    i = 0
    while i < len(input_grid):
        if input_grid[i] != 0 and input_grid[i] != 2:
            current_seq = []
            while i < len(input_grid) and input_grid[i] != 0 and input_grid[i] != 2:
                current_seq.append(input_grid[i])
                i += 1
            sequence = current_seq
        i += 1
    
    # Create output grid filled with zeros
    output = [0] * len(input_grid)
    
    # Place 2 in same position
    output[two_pos] = 2
    
    # Find where sequence is in input
    seq_start = -1
    for i in range(len(input_grid) - len(sequence) + 1):
        if all(input_grid[i+j] == sequence[j] for j in range(len(sequence))):
            seq_start = i
            break
    
    # If sequence is after 2, move it to the start
    if seq_start > two_pos:
        # Place sequence at start
        for i in range(len(sequence)):
            output[i] = sequence[i]
    else:
        # Keep sequence where it is
        for i in range(len(sequence)):
            output[seq_start + i] = sequence[i]
    
    return output

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0]
result = transform_grid(test_input)
print(' '.join(map(str, result)))