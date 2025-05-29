def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        current_num = 7
        next_num_pos = 0
        
        # First find all positions where numbers change
        changes = []
        for i in range(rows):
            if input_grid[i][col] not in [7, 9]:
                changes.append((i, input_grid[i][col]))
        
        # Fill the grid based on changes
        for i in range(rows):
            if changes and i >= changes[0][0]:
                current_num = changes[0][1]
                changes.pop(0)
            output[i][col] = current_num
            
            # Reset to 7 if we hit the next different number
            if changes and i + 1 == changes[0][0]:
                if current_num != changes[0][1]:
                    current_num = changes[0][1]
    
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