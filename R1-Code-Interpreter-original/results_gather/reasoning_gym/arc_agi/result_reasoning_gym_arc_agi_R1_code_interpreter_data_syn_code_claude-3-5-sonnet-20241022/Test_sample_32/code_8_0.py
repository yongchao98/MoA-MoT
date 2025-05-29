def transform_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Initialize output with first three columns
    output = [[rows[i][j] for j in range(3)] for i in range(height)]
    processed = [False] * height  # Track which rows have been processed
    
    # First pass: handle 4 patterns
    i = 0
    while i < height:
        if output[i][1] == 4:
            # Mark the next three rows (or remaining rows) as [4,4,4]
            for j in range(i, min(i+3, height)):
                output[j] = [4, 4, 4]
                processed[j] = True
            i += 3  # Skip to after the affected rows
        else:
            i += 1
    
    # Second pass: handle 9 patterns
    for i in range(height-1, -1, -1):  # Process from bottom up
        if not processed[i] and output[i][1] == 9:
            # Current row
            output[i] = [9, 9, 9]
            processed[i] = True
            
            # Row above (if exists and not processed)
            if i > 0 and not processed[i-1]:
                output[i-1] = [6, 9, 9]
                processed[i-1] = True
                
            # Two rows above (if exists and not processed)
            if i > 1 and not processed[i-2]:
                output[i-2] = [6, 6, 9]
                processed[i-2] = True
    
    # Keep original values for unprocessed rows
    for i in range(height):
        if not processed[i]:
            output[i] = [rows[i][0], rows[i][1], rows[i][2]]
    
    return '\n'.join(' '.join(str(x) for x in row) for row in output)

# Test with the input
test_input = """6 6 4 5 6 6 4 5 6 6 4 5 6 6 4
6 4 4 5 6 4 4 5 6 4 4 5 6 4 4
9 9 9 5 4 4 4 5 4 4 4 5 4 4 4
6 9 9 5 6 9 9 5 6 4 4 5 6 4 4
6 6 9 5 6 6 9 5 6 6 9 5 6 6 4
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
9 9 9 5 9 9 9 5 9 9 9 5 9 9 9
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
6 6 9 5 6 6 9 5 6 6 9 5 6 6 9
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
9 9 9 5 9 9 9 5 9 9 9 5 9 9 9
6 9 9 5 6 9 9 5 6 9 9 5 6 9 9
6 6 9 5 6 6 9 5 6 6 9 5 6 6 9"""

result = transform_grid(test_input)
print("Predicted output:")
print(result)