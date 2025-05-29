def transform_grid(input_grid):
    # Start with two 2s
    output_grid = [2, 2]
    
    # Count occurrences of 1s and 2s
    count_1 = input_grid.count(1)
    count_2 = input_grid.count(2)
    
    # Add 0s, 1s, and 2s to the output grid
    output_grid.extend([0] * (len(input_grid) - 2 - count_1 - count_2))
    output_grid.extend([1] * count_1)
    output_grid.extend([2] * (count_2 - 2))  # Subtract 2 because we already added two 2s at the start
    
    return output_grid

# Test input
test_input = [0, 1, 2, 1, 0, 0, 1, 0, 1, 0, 2]
test_output = transform_grid(test_input)
print(test_output)