def transform_grid(input_grid):
    # Find the first non-zero element
    first_non_zero_index = next((i for i, x in enumerate(input_grid) if x != 0), len(input_grid))
    
    # Create the output grid by shifting elements to the left
    shifted_elements = input_grid[first_non_zero_index:]
    output_grid = shifted_elements + [0] * (len(input_grid) - len(shifted_elements))
    
    return output_grid

# Test input
input_grid = [0, 0, 0, 0, 6, 3, 5, 1, 7, 9, 1, 2, 4, 9, 6, 3, 2, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)