def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        current_num = 7
        i = 0
        while i < rows:
            # If we find a new non-7 number
            if input_grid[i][col] not in [7, 9]:
                current_num = input_grid[i][col]
                start = i
                
                # Find where this section should end
                end = i + 1
                while end < rows:
                    # Stop if we find a different non-7 number
                    if input_grid[end][col] not in [7, 9] and input_grid[end][col] != current_num:
                        break
                    end += 1
                
                # Fill the section
                for j in range(start, end):
                    output[j][col] = current_num
                
                i = end
            else:
                i += 1
    
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