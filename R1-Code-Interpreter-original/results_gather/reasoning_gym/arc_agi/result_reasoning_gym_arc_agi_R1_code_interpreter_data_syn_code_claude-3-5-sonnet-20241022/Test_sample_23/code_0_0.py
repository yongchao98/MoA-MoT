def process_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [row[:] for row in input_grid]  # Create a copy
    
    # Find the middle section of zeros
    zero_section_start = None
    zero_section_end = None
    for i in range(rows):
        if all(x == 0 for x in input_grid[i]):
            if zero_section_start is None:
                zero_section_start = i
            zero_section_end = i
    
    # Process special numbers
    for col in range(cols):
        # Process upward from special numbers
        for row in range(zero_section_start-1, -1, -1):
            if input_grid[row][col] not in [0, 5]:
                # Create vertical path of zeros
                current_row = row
                while current_row < zero_section_start:
                    output[current_row][col] = 0
                    current_row += 1
        
        # Process downward from special numbers
        for row in range(zero_section_end+1, rows):
            if input_grid[row][col] not in [0, 5]:
                # Create vertical path of zeros
                current_row = row
                while current_row > zero_section_end:
                    output[current_row][col] = 0
                    current_row -= 1
    
    # Print the result
    for row in output:
        print(' '.join(map(str, row)))

# Test input
test_input = [
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,6],
    [5,5,5,5,5,5,6,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,6,5,5],
    [5,6,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,5,5,6,5,5,5,5,5,5,5],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,6,5,5,5,5,5,5,5,6,5,5],
    [5,6,5,5,5,5,5,5,5,5,5,5,5,5],
    [5,5,5,5,6,5,5,5,5,5,5,6,5,5],
    [5,5,5,5,5,5,5,5,5,5,5,5,5,5]
]

process_grid(test_input)