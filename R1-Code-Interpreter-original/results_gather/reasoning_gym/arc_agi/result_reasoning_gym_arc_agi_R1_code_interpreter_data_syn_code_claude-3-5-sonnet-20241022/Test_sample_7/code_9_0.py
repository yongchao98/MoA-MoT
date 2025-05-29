def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the unique number and first row where it appears
    unique_num = None
    first_pattern_row = 0
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != 4:
                unique_num = input_grid[i][j]
                first_pattern_row = i
                break
        if unique_num:
            break
    
    output = [[4] * cols for _ in range(rows)]
    
    # Calculate where the first group should start
    # If pattern appears in row 1, we start 2 rows before
    group_start = (first_pattern_row - 2) % 3
    
    # Process each row
    for i in range(rows):
        row_type = (i - group_start) % 3  # 0, 1, or 2
        
        # Fill from the start of each row
        if row_type == 0:  # First row of group: 9 4 4 pattern
            for j in range(0, cols, 3):
                if j < cols:
                    output[i][j] = unique_num
        elif row_type == 1:  # Second row of group: 9 4 4 pattern
            for j in range(0, cols, 3):
                if j < cols:
                    output[i][j] = unique_num
        else:  # Third row of group: 9 9 4 pattern
            for j in range(0, cols, 3):
                if j < cols:
                    output[i][j] = unique_num
                if j + 1 < cols:
                    output[i][j+1] = unique_num
    
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