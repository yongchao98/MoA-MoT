def reverse_segment(input_grid):
    # Find the first and last non-zero indices
    first_non_zero = None
    last_non_zero = None
    
    for i, value in enumerate(input_grid):
        if value != 0:
            if first_non_zero is None:
                first_non_zero = i
            last_non_zero = i
    
    # If no non-zero elements are found, return the input as is
    if first_non_zero is None or last_non_zero is None:
        return input_grid
    
    # Reverse the segment between the first and last non-zero indices
    output_grid = input_grid[:first_non_zero] + \
                  input_grid[first_non_zero:last_non_zero+1][::-1] + \
                  input_grid[last_non_zero+1:]
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Get the output grid
output_grid = reverse_segment(input_grid)

# Print the output in the required format
print("<<<" + " ".join(map(str, output_grid)) + ">>>")