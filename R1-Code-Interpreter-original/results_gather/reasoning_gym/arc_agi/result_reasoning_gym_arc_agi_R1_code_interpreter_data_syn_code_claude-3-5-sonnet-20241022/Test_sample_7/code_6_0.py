def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Find the unique number (non-4)
    unique_num = None
    for row in input_grid:
        for num in row:
            if num != 4:
                unique_num = num
                break
        if unique_num:
            break
    
    output = [[4] * cols for _ in range(rows)]
    
    # Process each row
    for i in range(rows):
        row_type = i % 3  # 0, 1, or 2
        
        # Fill from the start of each row
        if row_type == 2:  # Every third row
            for j in range(0, cols, 3):
                if j < cols:
                    output[i][j] = unique_num
                if j + 1 < cols:
                    output[i][j+1] = unique_num
        else:  # Other rows (0 and 1)
            for j in range(0, cols, 3):
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