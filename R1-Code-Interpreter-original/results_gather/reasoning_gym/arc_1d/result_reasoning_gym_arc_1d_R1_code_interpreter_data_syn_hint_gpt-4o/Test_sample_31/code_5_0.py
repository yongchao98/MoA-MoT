def transform_grid(input_grid):
    # Find the smallest non-zero number
    non_zero_numbers = [num for num in input_grid if num != 0]
    smallest_number = min(non_zero_numbers)
    
    # Create the output grid
    output_grid = [smallest_number] + [num for num in non_zero_numbers if num != smallest_number] + [0] * (len(input_grid) - len(non_zero_numbers))
    
    return output_grid

# Test input
input_grid = [4, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)