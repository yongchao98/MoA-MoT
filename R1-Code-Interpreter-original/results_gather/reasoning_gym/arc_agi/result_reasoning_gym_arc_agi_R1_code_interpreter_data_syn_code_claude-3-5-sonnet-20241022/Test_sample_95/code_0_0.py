def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column as it mostly remains unchanged
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        # Find sequences of non-7 numbers
        i = 0
        while i < rows:
            if input_grid[i][col] != 7:
                # Find the number to propagate
                num = input_grid[i][col]
                # Count sequence length
                seq_start = i
                while i < rows and (input_grid[i][col] != 7 or input_grid[i][col] == 9):
                    i += 1
                seq_end = i
                
                # Fill the sequence in output
                fill_start = max(0, seq_start)
                fill_length = seq_end - seq_start
                for j in range(fill_start, min(fill_start + fill_length, rows)):
                    output[j][col] = num
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