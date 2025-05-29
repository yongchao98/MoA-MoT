def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the unique number (non-4) and its leftmost position
    unique_num = None
    pattern_start = 0
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 4:
                unique_num = input_grid[i][j]
                pattern_start = j
                break
        if unique_num:
            break
    
    # Adjust pattern_start to be the starting position of the output pattern
    pattern_start = (pattern_start // 3) * 3
    
    output = [[4] * cols for _ in range(rows)]
    
    # Process each row
    for i in range(rows):
        row_type = i % 3  # 0, 1, or 2
        
        # Start filling from the pattern_start position
        if row_type == 2:  # Every third row (double pattern)
            for j in range(pattern_start, cols, 3):
                if j + 1 < cols:
                    output[i][j] = unique_num
                    output[i][j+1] = unique_num
        else:  # Other rows (single pattern)
            for j in range(pattern_start, cols, 3):
                if j < cols:
                    output[i][j] = unique_num
    
    # Print the output grid
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,9,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,9,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,9,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,4,4,9,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4],
    [4,4,4,4,4,4,4,4,9,9,4,4,4,4]
]

transform_grid(test_input)