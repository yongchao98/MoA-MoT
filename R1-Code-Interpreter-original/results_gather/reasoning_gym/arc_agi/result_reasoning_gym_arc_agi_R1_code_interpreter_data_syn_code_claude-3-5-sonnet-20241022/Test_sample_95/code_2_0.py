def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        # First, identify number groups
        groups = []  # (number, start_idx, end_idx)
        i = 0
        while i < rows:
            if input_grid[i][col] not in [7, 9]:
                num = input_grid[i][col]
                start = i
                # Find all connected instances of this number
                while i < rows and (input_grid[i][col] == num or input_grid[i][col] == 9):
                    i += 1
                end = i - 1
                groups.append((num, start, end))
            else:
                i += 1
        
        # Process each group
        for num, start, end in groups:
            # Fill from start to end
            for i in range(start, end + 1):
                output[i][col] = num
            
            # Extend upward if possible
            up = start - 1
            while up >= 0 and input_grid[up][col] == 7:
                output[up][col] = num
                up -= 1
                
            # Extend downward if possible
            down = end + 1
            while down < rows and input_grid[down][col] == 7:
                output[down][col] = num
                down += 1
    
    # Print the result
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [7, 0, 0],
    [7, 7, 9],
    [7, 7, 7],
    [7, 5, 9],
    [7, 9, 7],
    [7, 7, 7],
    [9, 6, 6],
    [7, 7, 9],
    [7, 7, 7],
    [7, 4, 9],
    [7, 4, 7]
]

process_grid(test_input)