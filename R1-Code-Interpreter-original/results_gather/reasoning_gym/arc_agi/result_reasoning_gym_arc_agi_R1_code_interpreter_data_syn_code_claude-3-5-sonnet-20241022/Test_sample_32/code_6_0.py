def transform_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Initialize output with first three columns
    output = [[rows[i][j] for j in range(3)] for i in range(height)]
    
    # First pass: identify 4 patterns
    for i in range(height):
        if output[i][1] == 4:  # Found a 4 in middle column
            # Affect current row and next two rows
            for j in range(i, min(i+3, height)):
                output[j] = [4, 4, 4]
    
    # Second pass: handle 9 patterns (but don't override 4 patterns)
    for i in range(height):
        if output[i][1] == 9 and output[i][0] != 4:  # Found a 9, not in 4 pattern
            if i > 0:  # Check row above
                output[i] = [6, 9, 9]
                if output[i-1][1] != 4:  # Don't override 4 pattern
                    output[i-1] = [6, 6, 9]
            else:
                output[i] = [6, 9, 9]
    
    # Keep original values where no patterns apply
    for i in range(height):
        if output[i][1] not in [4, 9] and output[i][1] != 4:
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