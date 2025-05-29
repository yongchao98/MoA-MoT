def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        # First pass: identify blocks and their numbers
        blocks = []  # List of (start, end, number)
        i = 0
        while i < rows:
            if input_grid[i][col] not in [7, 9]:
                num = input_grid[i][col]
                start = i
                # Look ahead for next different number
                i += 1
                while i < rows and input_grid[i][col] in [7, 9, num]:
                    i += 1
                blocks.append((start, i-1, num))
            else:
                i += 1
        
        # Second pass: fill blocks
        for block_idx, (start, end, num) in enumerate(blocks):
            # Determine next block start
            next_start = rows if block_idx == len(blocks)-1 else blocks[block_idx+1][0]
            # Fill until next block starts
            for i in range(start, next_start):
                output[i][col] = num
    
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