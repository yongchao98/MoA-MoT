def generate_output(input_str):
    # Convert input to list
    grid = [int(x) for x in input_str.split()]
    
    # Find position of 2 and sequence details
    two_pos = grid.index(2)
    sequence = [(i, x) for i, x in enumerate(grid) if x != 0 and x != 2]
    if not sequence:
        return grid
    
    sequence_num = sequence[0][1]
    sequence_count = len(sequence)
    
    # Create output grid
    output = [0] * 28
    
    # Place 2
    output[two_pos] = 2
    
    # Place sequence numbers right after 2
    for i in range(sequence_count):
        output[two_pos + 1 + i] = sequence_num
    
    # Ensure 3 trailing zeros
    output[-3:] = [0, 0, 0]
    
    return output

# Test with our input
test_input = "0 2 0 0 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0"
result = generate_output(test_input)
print(" ".join(str(x) for x in result))