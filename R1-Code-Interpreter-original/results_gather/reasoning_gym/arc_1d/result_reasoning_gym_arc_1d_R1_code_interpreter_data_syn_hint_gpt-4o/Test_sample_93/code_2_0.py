def transform_grid(input_grid):
    # Identify the sequence of non-zero numbers
    non_zero_sequence = [num for num in input_grid if num != 0]
    
    # Calculate the number of zeros to fill
    num_zeros = len(input_grid) - len(non_zero_sequence)
    
    # Construct the output grid
    output_grid = non_zero_sequence + [0] * num_zeros
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 2, 7, 0]
output = transform_grid(test_input)
print(output)