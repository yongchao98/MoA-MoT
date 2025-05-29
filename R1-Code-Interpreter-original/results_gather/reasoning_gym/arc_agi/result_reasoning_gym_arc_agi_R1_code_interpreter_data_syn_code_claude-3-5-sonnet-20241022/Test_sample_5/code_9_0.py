def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7 for _ in range(8)] for _ in range(rows)]
    out_row = 0
    
    # Process each 3-row segment
    for segment in range(0, rows-2, 3):
        # Get left side marker
        marker = None
        for r in range(segment, segment+3):
            for c in range(3):
                if input_grid[r][c] != 7 and input_grid[r][c] != 1:
                    marker = input_grid[r][c]
                    break
            if marker is not None:
                break
        
        # Set left side pattern
        if marker is not None:
            output[out_row][0] = marker
            output[out_row][1] = marker
        
        # Process right side patterns
        # First group (columns 7-8)
        if input_grid[segment][7] == 4 and input_grid[segment][8] != 4:
            value = input_grid[segment][8]
            output[out_row][2] = value
            output[out_row][3] = value
            
        # Second group (columns 10-11)
        if input_grid[segment][10] == 4 and input_grid[segment][11] != 4:
            value = input_grid[segment][11]
            output[out_row][4] = value
            output[out_row][5] = value
            
        # Third group (columns 13-14)
        if input_grid[segment][13] == 4 and input_grid[segment][14] != 4:
            value = input_grid[segment][14]
            output[out_row][6] = value
            output[out_row][7] = value
        
        out_row += 1
    
    # Print non-empty rows
    result = []
    for row in output:
        if any(x != 7 for x in row):
            result.append(' '.join(map(str, row)))
    print('\n'.join(result))

# Input grid
input_text = """7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 6 7 7 7 7 7 4 4 7 4 4 7 4 4 7 4 4 7 7
7 7 1 7 7 7 7 4 2 7 4 4 7 4 9 7 4 9 7 7
7 7 1 1 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 7 7 7 7 7 4 4 7 4 4 7 4 4 7 4 4 7 7
7 2 7 7 7 7 7 4 2 7 4 6 7 4 6 7 4 4 7 7
7 7 1 1 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 1 7 7 7 7 4 4 7 4 4 7 4 4 7 4 4 7 7
7 7 7 7 7 7 7 4 2 7 4 6 7 4 4 7 4 4 7 7
7 0 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 1 1 7 7 7 4 4 7 4 4 7 4 4 7 4 4 7 7
7 7 1 1 7 7 7 4 0 7 4 4 7 4 4 7 4 0 7 7
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 9 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 7 1 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 1 1 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7
7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7"""

input_grid = [list(map(int, line.split())) for line in input_text.split('\n')]
process_grid(input_grid)