def transform_grid(input_grid):
    # Find the longest contiguous sequence of non-zero numbers
    max_sequence = []
    current_sequence = []
    
    for num in input_grid:
        if num != 0:
            current_sequence.append(num)
        else:
            if len(current_sequence) > len(max_sequence):
                max_sequence = current_sequence
            current_sequence = []
    
    # Check at the end of the loop
    if len(current_sequence) > len(max_sequence):
        max_sequence = current_sequence
    
    # Determine the number of leading zeros
    leading_zeros = 0
    for num in input_grid:
        if num == 0:
            leading_zeros += 1
        else:
            break
    
    # Create the output grid
    output_grid = [0] * leading_zeros + max_sequence + [0] * (len(input_grid) - leading_zeros - len(max_sequence))
    
    return output_grid

# Test input
test_input = [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
output = transform_grid(test_input)
print(output)