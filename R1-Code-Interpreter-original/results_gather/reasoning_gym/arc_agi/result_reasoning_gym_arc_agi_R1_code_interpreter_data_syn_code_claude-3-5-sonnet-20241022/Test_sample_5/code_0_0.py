def process_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Initialize output grid (8 columns as seen in examples)
    output = [[7 for _ in range(8)] for _ in range(11)]
    
    # Process left side patterns
    for i in range(0, rows-2, 3):
        if i + 2 >= rows:
            break
            
        # Get the left-side marker (first non-7 in the rows)
        marker = None
        marker_row = None
        for r in range(i, i+3):
            for c in range(3):
                if input_grid[r][c] != 7:
                    marker = input_grid[r][c]
                    marker_row = r
                    break
            if marker is not None:
                break
        
        if marker is None:
            continue
            
        # Count pattern of 1's
        ones_pattern = []
        for r in range(i, min(i+3, rows)):
            count = 0
            for c in range(3):
                if input_grid[r][c] == 1:
                    count += 1
            if count > 0:
                ones_pattern.append(count)
        
        out_row = i // 3
        if len(ones_pattern) > 0:
            # Fill left side based on marker and ones pattern
            if len(ones_pattern) >= 2:
                output[out_row][0] = marker
                output[out_row][1] = marker
                if len(ones_pattern) > 2:
                    output[out_row+1][0] = marker
                    output[out_row+1][1] = 7
            else:
                output[out_row][0] = marker
                output[out_row][1] = 7
    
    # Process right side patterns
    for i in range(0, rows-2, 3):
        if i + 2 >= rows:
            break
            
        # Look for patterns in the 4's sections
        for section in range(3):  # There are 3 sections of 4's
            base_col = 7 + section * 3
            if base_col + 1 >= cols:
                break
                
            if input_grid[i][base_col] == 4:
                marker = input_grid[i][base_col + 1]
                out_col = 2 + section * 2
                
                # Check pattern in this section
                if marker != 4:
                    output[i//3][out_col] = marker
                    output[i//3][out_col+1] = marker
                    
                    # Check next row
                    next_marker = input_grid[i+1][base_col + 1]
                    if next_marker != marker and next_marker != 4:
                        output[i//3+1][out_col] = marker
                        output[i//3+1][out_col+1] = next_marker
    
    # Print the output grid
    result = []
    for row in output:
        if any(x != 7 for x in row):  # Only print rows that have non-7 values
            result.append(' '.join(map(str, row)))
    print('\n'.join(result))

# Parse the input grid
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