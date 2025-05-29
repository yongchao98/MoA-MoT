def analyze_pattern(input_grid, output_grid):
    # Convert string grids to 2D arrays
    def parse_grid(grid):
        return [[int(x) for x in row.split()] for row in grid.strip().split('\n')]
    
    def find_patterns(grid):
        rows, cols = len(grid), len(grid[0])
        patterns = []
        
        # Check horizontal patterns (3 consecutive 9s)
        for i in range(rows):
            for j in range(cols-2):
                if grid[i][j:j+3] == [9,9,9]:
                    patterns.append(('H', i, j))
        
        # Check vertical patterns (3 consecutive 9s)
        for i in range(rows-2):
            for j in range(cols):
                if [grid[i][j], grid[i+1][j], grid[i+2][j]] == [9,9,9]:
                    patterns.append(('V', i, j))
        
        return patterns

    # Test input
    test_input = """9 9 9 8 8 8 8 8 8
8 9 8 8 8 8 8 9 8
8 8 8 9 8 8 9 9 9
8 8 9 9 9 8 8 8 8
8 8 8 8 8 9 8 8 8
9 9 9 8 9 9 9 8 8
8 9 8 8 8 8 8 8 8
8 8 8 8 8 8 9 9 9
8 8 9 9 9 8 8 9 8
8 8 8 9 8 8 8 8 8"""
    
    input_arr = [[int(x) for x in row.split()] for row in test_input.split('\n')]
    patterns = find_patterns(input_arr)
    
    # Generate output based on patterns
    rows = len(input_arr) // 3
    cols = len(input_arr[0]) // 2
    output = [[8]*cols for _ in range(rows)]
    
    # For each pattern found, mark corresponding position in output
    for pattern_type, row, col in patterns:
        out_row = row // 3
        out_col = col // 2
        if pattern_type == 'H':
            output[out_row][out_col] = 1
            if out_col + 1 < cols:
                output[out_row][out_col + 1] = 3
        elif pattern_type == 'V':
            output[out_row][out_col] = 1
            if out_row + 1 < rows:
                output[out_row + 1][out_col] = 3
    
    # Print the output
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = """9 9 9 8 8 8 8 8 8
8 9 8 8 8 8 8 9 8
8 8 8 9 8 8 9 9 9
8 8 9 9 9 8 8 8 8
8 8 8 8 8 9 8 8 8
9 9 9 8 9 9 9 8 8
8 9 8 8 8 8 8 8 8
8 8 8 8 8 8 9 9 9
8 8 9 9 9 8 8 9 8
8 8 8 9 8 8 8 8 8"""

analyze_pattern(test_input, None)