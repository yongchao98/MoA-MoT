def find_output_pattern(input_grid):
    # Find all non-6 numbers and their positions
    numbers = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 6:
                numbers.append((i, j, input_grid[i][j]))
    
    # Sort numbers by position
    # First by column (right to left), then by row (top to bottom)
    numbers.sort(key=lambda x: (-x[1], x[0]))
    
    # Initialize output grid
    output = [[0 for _ in range(3)] for _ in range(3)]
    
    # Fill output grid based on the pattern observed
    # First row: rightmost numbers
    # Second row: middle numbers
    # Third row: leftmost numbers
    for i, (row, col, val) in enumerate(numbers):
        if col >= 7:  # rightmost section
            output_col = 2
            output_row = 0
        elif col >= 3:  # middle section
            output_col = 1
            output_row = 1
        else:  # leftmost section
            output_col = 0
            output_row = 2
            
        # Adjust position based on row in input
        if row >= 7:
            output_col = min(2, output_col + 1)
        elif row >= 3:
            output_col = output_col
        else:
            output_col = max(0, output_col - 1)
            
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