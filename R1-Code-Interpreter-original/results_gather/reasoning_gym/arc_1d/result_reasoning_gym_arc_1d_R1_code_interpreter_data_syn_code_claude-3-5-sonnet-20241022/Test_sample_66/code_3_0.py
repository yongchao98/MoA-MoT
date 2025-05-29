def generate_final_output(input_str):
    # Convert input to list
    grid = [int(x) for x in input_str.split()]
    output = [0] * 28  # Initialize with 28 zeros
    
    # Find position of 2 and sequence details
    two_pos = grid.index(2)
    sequence = [(i, x) for i, x in enumerate(grid) if x != 0 and x != 2]
    if not sequence:
        return grid
    
    sequence_num = sequence[0][1]
    sequence_count = len(sequence)
    
    # Place 2
    output[two_pos] = 2
    
    # Calculate start position for sequence to center it
    available_space = 28 - 3 - (two_pos + 1)  # Total - trailing zeros - (2's position + 1)
    sequence_start = two_pos + 1
    
    # Place sequence numbers
    for i in range(sequence_count):
        output[sequence_start + i] = sequence_num
    
    # Ensure last 3 positions are zeros
    output[-3:] = [0, 0, 0]
    
    return output

# Test with our input
test_input = "0 2 0 0 0 0 0 0 0 0 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0"
result = generate_final_output(test_input)
print(" ".join(str(x) for x in result))