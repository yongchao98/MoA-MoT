def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7]*3 for _ in range(rows)]
    
    # Copy third column
    for i in range(rows):
        output[i][2] = input_grid[i][2]
    
    # Process first two columns
    for col in range(2):
        # First find all non-7, non-9 numbers and their positions
        number_positions = []
        for i in range(rows):
            if input_grid[i][col] not in [7, 9]:
                number_positions.append((i, input_grid[i][col]))
        
        # Process each section
        for i in range(len(number_positions)):
            start_pos, num = number_positions[i]
            # Find end position (either next different number or end of grid)
            end_pos = rows
            if i + 1 < len(number_positions):
                next_pos, next_num = number_positions[i + 1]
                if next_num != num:
                    end_pos = next_pos
            
            # Fill the section
            for j in range(start_pos, end_pos):
                output[j][col] = num
    
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