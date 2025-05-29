def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        i = 0
        while i < rows:
            if input_grid[i][col] not in [7, 9]:
                current_num = input_grid[i][col]
                # Find the end of this block
                start = i
                while i < rows:
                    if input_grid[i][col] not in [7, 9, current_num]:
                        break
                    i += 1
                # Fill the block
                for j in range(start, i):
                    output[j][col] = current_num
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