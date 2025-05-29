def transform_grid(input_grid):
    output_grid = input_grid[:]
    
    # Extend leading non-zero numbers
    for i in range(len(input_grid)):
        if input_grid[i] != 0:
            leading_number = input_grid[i]
            break
    for j in range(i):
        output_grid[j] = leading_number
    
    # Extend trailing non-zero numbers
    for i in range(len(input_grid) - 1, -1, -1):
        if input_grid[i] != 0:
            trailing_number = input_grid[i]
            break
    for j in range(i, len(input_grid)):
        output_grid[j] = trailing_number
    
    # Extend middle '1's to a length of four
    middle_start = input_grid.index(1)
    middle_end = middle_start + input_grid[middle_start:].count(1)
    for j in range(middle_start, middle_start + 4):
        if j < len(output_grid):
            output_grid[j] = 1
    
    return output_grid

# Test input
test_input = [0, 0, 0, 0, 7, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
test_output = transform_grid(test_input)
print(test_output)