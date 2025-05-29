def transform_grid(input_grid):
    # Find the majority number
    majority_number = max(set(input_grid), key=input_grid.count)
    
    # Find the first non-majority number
    for num in input_grid:
        if num != majority_number and num != 0:
            non_majority_number = num
            break
    
    # Find the position of the first majority number
    first_majority_index = input_grid.index(majority_number)
    
    # Create the output grid
    output_grid = input_grid[:]
    output_grid[first_majority_index] = non_majority_number
    
    # Shift the rest of the sequence to the right
    for i in range(first_majority_index + 1, len(input_grid)):
        if input_grid[i] == non_majority_number:
            output_grid[i] = majority_number
        else:
            output_grid[i] = input_grid[i - 1]
    
    # Ensure the majority number is placed after the moved number
    output_grid[first_majority_index + 1] = majority_number
    
    return output_grid

# Test input
test_input = [4, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
output = transform_grid(test_input)
print(output)