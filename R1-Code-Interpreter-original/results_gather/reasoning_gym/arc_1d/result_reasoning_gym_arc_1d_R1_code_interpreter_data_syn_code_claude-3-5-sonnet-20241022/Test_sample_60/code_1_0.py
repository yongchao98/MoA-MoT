def transform_grid(input_grid):
    output = [0] * len(input_grid)
    n = len(input_grid)
    
    # Find all sequences first
    sequences = []
    i = 0
    while i < n:
        if input_grid[i] != 0:
            start = i
            value = input_grid[i]
            while i < n and input_grid[i] != 0:
                i += 1
            sequences.append((start, value))
        else:
            i += 1
    
    # Apply patterns at fixed intervals
    pattern_start = 0
    seq_index = 0
    
    while pattern_start < n and seq_index < len(sequences):
        # Fill pattern
        for i in range(pattern_start, min(n, pattern_start + 5)):
            output[i] = sequences[seq_index][1]
        
        # Move to next pattern position
        pattern_start += 5
        seq_index += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Print result
print(transform_grid(test_input))