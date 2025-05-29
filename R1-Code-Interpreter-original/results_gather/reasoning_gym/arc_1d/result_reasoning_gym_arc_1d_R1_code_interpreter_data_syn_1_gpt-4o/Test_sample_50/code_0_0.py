def transform_grid(input_grid):
    # Initialize the output grid with zeros
    output_grid = [0] * 25
    
    # Extract the non-zero elements from the input grid
    non_zero_elements = [x for x in input_grid if x != 0]
    
    # Place the non-zero elements at the start of the output grid
    for i in range(len(non_zero_elements)):
        output_grid[i] = non_zero_elements[i]
    
    # Return the transformed output grid
    return output_grid

# Test input grid
test_input = [2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

# Apply the transformation
test_output = transform_grid(test_input)

# Print the output grid
print(test_output)