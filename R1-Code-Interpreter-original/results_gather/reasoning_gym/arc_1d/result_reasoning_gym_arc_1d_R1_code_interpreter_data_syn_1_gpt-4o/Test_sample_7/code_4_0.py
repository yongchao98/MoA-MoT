def transform_grid(input_grid):
    # Find the non-zero elements
    non_zero_elements = [x for x in input_grid if x != 0]
    
    # Calculate the number of zeros to add on each side
    total_length = len(input_grid)
    non_zero_length = len(non_zero_elements)
    zeros_to_add = total_length - non_zero_length
    
    # Calculate left and right zeros
    left_zeros = zeros_to_add // 2
    right_zeros = zeros_to_add - left_zeros
    
    # Create the output grid
    output_grid = [0] * left_zeros + non_zero_elements + [0] * right_zeros
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)