def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7 for _ in range(8)] for _ in range(rows)]
    out_row = 0
    
    # Process each 3-row segment
    for segment in range(0, rows-2, 3):
        # Process left side pattern
        marker = None
        for r in range(segment, segment+3):
            for c in range(3):
                if input_grid[r][c] != 7 and input_grid[r][c] != 1:
                    marker = input_grid[r][c]
                    break
            if marker is not None:
                break
        
        if marker is not None:
            output[out_row][0] = marker
            output[out_row][1] = marker
        
        # Process right side patterns (three groups)
        for group in range(3):
            base_col = 7 + group * 3
            out_col = 2 + group * 2
            
            if base_col + 1 >= len(input_grid[0]):
                continue
            
            # First row pattern
            if input_grid[segment][base_col] == 4:
                value = input_grid[segment][base_col + 1]
                if value != 4:
                    output[out_row][out_col] = value
                    output[out_row][out_col + 1] = value
            
            # Second row pattern
            if input_grid[segment+1][base_col] == 4:
                value = input_grid[segment+1][base_col + 1]
                if value != 4:
                    output[out_row][out_col] = value
                    output[out_row][out_col + 1] = value
            
            # Third row pattern (if exists)
            if segment+2 < rows and input_grid[segment+2][base_col] == 4:
                value = input_grid[segment+2][base_col + 1]
                if value != 4:
                    output[out_row][out_col] = value
                    output[out_row][out_col + 1] = value
        
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