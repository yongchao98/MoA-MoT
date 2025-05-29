def find_output_pattern(input_grid):
    # Find all non-6 numbers and their positions
    numbers = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 6:
                numbers.append((i, j, input_grid[i][j]))
    
    # Sort numbers by column first (right to left), then by row (top to bottom)
    numbers.sort(key=lambda x: (-x[1], x[0]))
    
    # Initialize output grid
    output = [[0 for _ in range(3)] for _ in range(3)]
    
    # Map first three numbers to first row (right to left)
    row_idx = 0
    col_idx = 2
    for i, (_, _, val) in enumerate(numbers):
        if i < 3:
            output[0][col_idx] = val
            col_idx -= 1
        elif i < 6:
            output[1][col_idx+1] = val
            col_idx -= 1
        else:
            output[2][col_idx+1] = val
            col_idx -= 1
        if col_idx < 0:
            col_idx = 2
    
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