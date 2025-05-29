def transform_grid(input_grid):
    # Extract non-zero elements
    non_zero_elements = [num for num in input_grid if num != 0]
    # Calculate the number of zeros to append
    num_zeros = len(input_grid) - len(non_zero_elements)
    # Create the output grid
    output_grid = non_zero_elements + [0] * num_zeros
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 1, 5, 0, 0, 0]
output_grid = transform_grid(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")