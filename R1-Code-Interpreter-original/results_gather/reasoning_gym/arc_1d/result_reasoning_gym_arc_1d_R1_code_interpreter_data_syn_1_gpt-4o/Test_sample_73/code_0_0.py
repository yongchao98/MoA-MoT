def transform_grid(input_grid):
    # Find the start and end of the sequence of non-zero numbers
    start = None
    end = None
    for i in range(len(input_grid)):
        if input_grid[i] != 0:
            if start is None:
                start = i
            end = i
    
    # If a sequence is found, shift it one position to the left
    if start is not None and end is not None:
        output_grid = input_grid[:]
        for i in range(start, end):
            output_grid[i] = input_grid[i + 1]
        output_grid[end] = 0
        return output_grid
    else:
        return input_grid

# Test input
test_input = [0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 2, 0]
output = transform_grid(test_input)
print(output)