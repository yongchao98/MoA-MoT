def find_output_pattern(input_grid):
    # Find all non-6 numbers and their positions
    numbers = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 6:
                numbers.append((i, j, input_grid[i][j]))
    
    # Initialize output grid
    output = [[6 for _ in range(3)] for _ in range(3)]
    
    # Sort numbers by their position in a specific way
    # First by row (top to bottom), then by column (right to left)
    numbers.sort(key=lambda x: (x[0], -x[1]))
    
    # Map the numbers to output grid
    for i, (row, col, val) in enumerate(numbers):
        # Determine output position based on input position
        if row < 4:  # top third
            output_row = 0
        elif row < 7:  # middle third
            output_row = 1
        else:  # bottom third
            output_row = 2
            
        if col < 4:  # left third
            output_col = 0
        elif col < 7:  # middle third
            output_col = 1
        else:  # right third
            output_col = 2
            
        output[output_row][output_col] = val
    
    # Print output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [6, 6, 0, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 4, 6, 6, 6],
    [3, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 1, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 0, 6],
    [6, 5, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 3, 6],
    [6, 6, 6, 6, 6, 2, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 9, 6, 6]
]

find_output_pattern(test_input)