def process_grid(input_grid):
    rows = len(input_grid)
    output = [[7 for _ in range(8)] for _ in range(15)]  # Initialize larger than needed
    out_row = 0
    
    # Process 3-row segments
    for segment in range(0, rows-2, 3):
        segment_rows = input_grid[segment:segment+3]
        
        # Find left side marker and 1's pattern
        marker = None
        ones_count = [0, 0, 0]
        for i, row in enumerate(segment_rows):
            for j in range(3):
                if row[j] != 7:
                    if row[j] == 1:
                        ones_count[i] += 1
                    elif marker is None:
                        marker = row[j]
        
        # Process left side pattern
        if marker is not None:
            output[out_row][0] = marker
            output[out_row][1] = marker
            if ones_count[1] > 0 or ones_count[2] > 0:
                output[out_row+1][0] = marker
                if ones_count[2] > 0:
                    output[out_row+1][1] = marker
        
        # Process right side patterns (4's groups)
        for group in range(3):
            base_col = 7 + group * 3
            out_col = 2 + group * 2
            
            if base_col + 1 >= len(input_grid[0]):
                break
                
            # Check for 4-X pattern
            if segment_rows[0][base_col] == 4:
                value = segment_rows[0][base_col + 1]
                if value != 4:
                    output[out_row][out_col] = value
                    output[out_row][out_col + 1] = value
                    
                    # Check next row
                    if segment_rows[1][base_col] == 4:
                        next_value = segment_rows[1][base_col + 1]
                        if next_value != 4 and next_value != value:
                            output[out_row+1][out_col] = value
                            output[out_row+1][out_col + 1] = next_value
        
        out_row += 3
    
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