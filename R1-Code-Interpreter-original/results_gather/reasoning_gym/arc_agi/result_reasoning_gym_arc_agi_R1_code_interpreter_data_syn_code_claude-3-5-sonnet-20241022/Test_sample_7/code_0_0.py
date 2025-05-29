def analyze_and_transform(input_grid):
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
    
    # Create output grid
    output = [[4] * cols for _ in range(rows)]
    
    # For each row in the output
    for i in range(rows):
        # Find pattern in input row
        pattern_positions = []
        for j in range(cols):
            if input_grid[i][j] == unique_num:
                pattern_positions.append(j)
        
        # If this row contains the unique number
        if pattern_positions:
            # Fill the entire row with pattern
            for j in range(cols):
                # Every third position gets a 4, others get unique_num
                if (j % 3) == 1:
                    output[i][j] = 4
                else:
                    output[i][j] = unique_num
        
        # If row above had pattern
        elif i > 0 and any(input_grid[i-1][j] == unique_num for j in range(cols)):
            # Fill with alternating pattern
            for j in range(cols):
                if j % 3 == 1:
                    output[i][j] = 4
                else:
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

analyze_and_transform(test_input)