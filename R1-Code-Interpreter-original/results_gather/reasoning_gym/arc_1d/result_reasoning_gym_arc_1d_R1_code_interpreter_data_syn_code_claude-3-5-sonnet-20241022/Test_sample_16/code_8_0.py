def transform_grid(input_grid):
    # Find sequences in the input
    sequences = []
    current_seq = [input_grid[0]]
    
    # Group numbers into sequences
    for i in range(1, len(input_grid)):
        if (input_grid[i] == 0 and input_grid[i-1] == 0) or \
           (input_grid[i] != 0 and input_grid[i-1] != 0):
            current_seq.append(input_grid[i])
        else:
            sequences.append(current_seq)
            current_seq = [input_grid[i]]
    sequences.append(current_seq)
    
    # Transform based on sequences
    output = []
    
    # If starts with exactly two zeros
    if len(sequences[0]) == 2 and sequences[0][0] == 0 and sequences[0][1] == 0:
        output = [0, 0, 0]  # Add three zeros
        # Add remaining sequences except last element
        for seq in sequences[1:]:
            if seq[-1] == input_grid[-1]:  # If it's the last sequence
                output.extend(seq[:-1])
            else:
                output.extend(seq)
    else:
        output = [0] + input_grid[:-1]
    
    return output

# Test all examples
example1_in = [2, 1, 0, 0, 9, 4, 0, 0, 3, 0, 0, 0, 0, 0]
example1_out = [0, 2, 1, 0, 0, 9, 4, 0, 0, 3, 0, 0, 0, 0]

example2_in = [0, 0, 7, 7, 3, 0, 0, 1, 4, 0, 6, 9, 0, 3]
example2_out = [0, 0, 0, 7, 7, 3, 0, 0, 1, 4, 0, 6, 9, 3]

example3_in = [5, 7, 6, 0, 1, 0, 6, 6, 0, 0, 0, 5, 2, 0]
example3_out = [0, 5, 7, 6, 0, 1, 0, 6, 6, 0, 0, 0, 5, 2]

# Verify examples
print("Example 1:", transform_grid(example1_in) == example1_out)
print("Example 2:", transform_grid(example2_in) == example2_out)
print("Example 3:", transform_grid(example3_in) == example3_out)

# Transform test input
test_input = [0, 8, 0, 0, 0, 0, 0, 0, 7, 3, 0, 0, 4, 6]
result = transform_grid(test_input)
print("\nTest input:", test_input)
print("Result:", result)