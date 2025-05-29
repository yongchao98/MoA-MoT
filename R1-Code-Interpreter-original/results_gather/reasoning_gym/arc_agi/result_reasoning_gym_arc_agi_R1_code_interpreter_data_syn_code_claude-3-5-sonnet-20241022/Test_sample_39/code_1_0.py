def transform_grid(input_grid):
    # Convert input strings to 2D array
    grid = [row.split() for row in input_grid]
    output = [row[:] for row in grid]
    
    # Process every other row starting from row 0
    for i in range(0, len(grid)-2, 2):
        # Find the 0 in the current group of three rows
        row1 = grid[i]
        row2 = grid[i+1]
        row3 = grid[i+2] if i+2 < len(grid) else None
        
        # If first row has a 0, check pattern
        if '0' in row1:
            pos = row1.index('0')
            if row2[pos] == '9' and (row3 is None or row3[pos] == '9'):
                output[i][pos] = '5'
        
        # If middle row has a 0, check pattern
        if row2.count('0') == 1:
            pos = row2.index('0')
            if row1[pos] == '9' and (row3 is None or row3[pos] == '9'):
                output[i+1][pos] = '5'
        
        # If last row exists and has a 0, check pattern
        if row3 is not None and '0' in row3:
            pos = row3.index('0')
            if row1[pos] == '9' and row2[pos] == '9':
                output[i+2][pos] = '5'
    
    # Convert back to string format
    return [' '.join(row) for row in output]

# Test input
test_input = [
    "0 9 9", "9 0 9", "9 9 0", "0 9 9", "9 9 0", "9 0 9",
    "0 9 9", "9 9 0", "9 0 9", "9 9 0", "0 9 9", "9 0 9",
    "0 9 9", "9 9 0", "9 0 9", "0 9 9", "9 0 9"
]

result = transform_grid(test_input)
for row in result:
    print(row)