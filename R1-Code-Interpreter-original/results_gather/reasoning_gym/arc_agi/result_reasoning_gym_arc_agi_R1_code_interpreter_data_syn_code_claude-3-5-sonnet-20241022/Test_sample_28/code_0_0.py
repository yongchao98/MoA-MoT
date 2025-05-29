def analyze_grid(input_grid):
    # Find positions of non-6 numbers
    positions = []
    for i in range(10):
        for j in range(10):
            if input_grid[i][j] != 6:
                positions.append((i, j, input_grid[i][j]))
    
    # Sort positions by row and column to analyze pattern
    positions.sort()
    
    # Create output grid
    output = [[0 for _ in range(3)] for _ in range(3)]
    
    # Analyze positions and map to output
    for pos in positions:
        row, col, value = pos
        # Divide the 10x10 grid into 3x3 regions
        output_row = row // 4
        output_col = col // 4
        if output_row >= 3:
            output_row = 2
        if output_col >= 3:
            output_col = 2
        output[output_row][output_col] = value
    
    # Print the output grid
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

analyze_grid(test_input)