def find_sequence(grid):
    # Find non-zero sequence (excluding 2)
    sequence = []
    current_seq = []
    for num in grid:
        if num != 0 and num != 2:
            current_seq.append(num)
        elif len(current_seq) > 0:
            sequence = current_seq
            current_seq = []
    if len(current_seq) > 0:
        sequence = current_seq
    return sequence

def find_two_position(grid):
    return grid.index(2)

def transform_grid(input_grid):
    # Find position of 2
    two_pos = find_two_position(input_grid)
    
    # Find the sequence
    sequence = find_sequence(input_grid)
    
    # Create output grid filled with zeros
    output = [0] * len(input_grid)
    
    # Place 2 in same position
    output[two_pos] = 2
    
    # Find where sequence is in input (before or after 2)
    seq_start = -1
    for i in range(len(input_grid) - len(sequence) + 1):
        if all(input_grid[i+j] == sequence[j] for j in range(len(sequence))):
            seq_start = i
            break
    
    # If sequence is after 2, move it before 2
    if seq_start > two_pos:
        # Place sequence before 2
        for i in range(len(sequence)):
            output[two_pos - len(sequence) + i] = sequence[i]
    else:
        # Keep sequence where it is
        for i in range(len(sequence)):
            output[seq_start + i] = sequence[i]
    
    return output

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0]
result = transform_grid(test_input)
print(' '.join(map(str, result)))