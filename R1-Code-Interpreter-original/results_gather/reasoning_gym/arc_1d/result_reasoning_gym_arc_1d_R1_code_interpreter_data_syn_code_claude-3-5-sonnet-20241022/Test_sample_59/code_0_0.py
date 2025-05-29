def analyze_and_transform(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Find the non-zero sequence
    non_zero_start = -1
    non_zero_end = -1
    non_zero_value = 0
    
    for i in range(len(grid)):
        if grid[i] != 0 and non_zero_start == -1:
            non_zero_start = i
            non_zero_value = grid[i]
        elif grid[i] == 0 and non_zero_start != -1 and non_zero_end == -1:
            non_zero_end = i - 1
            break
    
    if non_zero_end == -1:
        non_zero_end = len(grid) - 1
        
    # Calculate sequence length
    sequence_length = non_zero_end - non_zero_start + 1
    
    # Create output grid with left shift of 3
    output = [0] * len(grid)
    new_start = non_zero_start - 3
    
    # Fill in the shifted sequence
    for i in range(sequence_length):
        output[new_start + i] = non_zero_value
    
    # Convert to string format
    return ' '.join(map(str, output))

# Test input
test_input = "0 0 0 0 0 0 0 0 0 0 8 8 8 8 8 8 8 0 0 0 0 0 0 0 0 0 0 0"
result = analyze_and_transform(test_input)
print(result)