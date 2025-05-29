def transform_grid(input_grid):
    # Copy the first five elements
    output_grid = input_grid[:5]
    
    # Initialize the rest of the grid with zeros
    rest_of_grid = [0] * (len(input_grid) - 5)
    
    # Fill the rest of the grid following the pattern observed
    # The pattern seems to involve filling the zeros with the non-zero number
    # in a specific way, creating a diagonal or a specific pattern.
    non_zero_number = input_grid[0]
    
    # Fill the pattern
    rest_of_grid[1] = non_zero_number
    rest_of_grid[2:5] = [non_zero_number] * 3
    rest_of_grid[6:9] = [non_zero_number] * 3
    rest_of_grid[10] = non_zero_number
    
    # Combine the two parts
    output_grid.extend(rest_of_grid)
    
    return output_grid

# Test input
test_input = [1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)