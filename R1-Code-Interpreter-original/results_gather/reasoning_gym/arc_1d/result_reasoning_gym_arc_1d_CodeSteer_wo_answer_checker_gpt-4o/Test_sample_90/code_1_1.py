def transform_grid(input_grid):
    # Create a copy of the input grid to avoid modifying the original
    output_grid = input_grid[:]
    
    # Apply the transformation rule
    if len(output_grid) > 1:
        output_grid[1] = 0  # Set the second element to 0
    
    return output_grid

# Test input
test_input = [9, 9, 9, 0, 9, 0, 9, 9, 0, 9, 9, 9, 9, 0, 0, 0]

# Transform the test input using the function
test_output = transform_grid(test_input)

# Print the output in the required format
print("<<<" + " ".join(map(str, test_output)) + ">>>")