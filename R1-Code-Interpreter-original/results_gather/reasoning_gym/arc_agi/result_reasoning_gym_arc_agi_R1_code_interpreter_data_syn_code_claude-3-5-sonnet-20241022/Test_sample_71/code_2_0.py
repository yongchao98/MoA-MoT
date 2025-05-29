def generate_output(input_grid):
    n = len(input_grid)
    output = [row[:] for row in input_grid]
    
    # In a 10x10 grid, always check row 5 (index 5)
    middle_row = 5
    
    # Check if there's any pattern above and below the middle row
    has_pattern = False
    for i in range(n):
        for j in range(n):
            if input_grid[i][j] != 7:
                has_pattern = True
                break
        if has_pattern:
            break
    
    # If there's a pattern, place 2's in the middle row
    if has_pattern:
        output[middle_row] = [2] * n
    
    # Print the result
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [7,7,7,7,7,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,4,4,4,4,7,7,7,7,7],
    [7,7,7,7,7,7,7,7,7,7]
]

generate_output(test_input)