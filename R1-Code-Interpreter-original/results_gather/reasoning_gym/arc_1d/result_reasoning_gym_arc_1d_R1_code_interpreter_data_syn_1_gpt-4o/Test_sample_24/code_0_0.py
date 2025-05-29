def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), None)
    if first_non_zero_index is None:
        return input_grid  # No non-zero elements, return as is

    # Find the last non-zero element
    last_non_zero_index = len(input_grid) - next((i for i, x in enumerate(reversed(input_grid)) if x != 0), None) - 1

    # Extract the sequence of non-zero elements
    non_zero_sequence = input_grid[first_non_zero_index:last_non_zero_index + 1]

    # Calculate the number of leading zeros in the input
    leading_zeros = first_non_zero_index

    # Calculate the new number of leading zeros
    new_leading_zeros = leading_zeros + 2

    # Create the new output grid
    output_grid = [0] * new_leading_zeros + non_zero_sequence + [0] * (len(input_grid) - new_leading_zeros - len(non_zero_sequence))
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)