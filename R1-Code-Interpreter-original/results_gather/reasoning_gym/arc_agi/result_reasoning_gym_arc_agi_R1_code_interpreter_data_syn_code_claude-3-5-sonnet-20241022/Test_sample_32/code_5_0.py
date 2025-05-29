def transform_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Initialize output with first three columns
    output = [[rows[i][j] for j in range(3)] for i in range(height)]
    
    # First pass: handle 4 patterns (they take precedence)
    for i in range(height-1, -1, -1):
        if output[i][1] == 4:
            # If we find a 4 in middle column, affect current and two rows above
            for j in range(max(0, i-2), i+1):
                output[j] = [4, 4, 4]
    
    # Second pass: handle 9 patterns (but don't override 4 patterns)
    for i in range(height-1, -1, -1):
        if output[i][1] == 9 and output[i][0] != 4:  # Don't override 4 patterns
            if i < height-1 and output[i+1][1] == 9:  # Check if part of a sequence
                output[i] = [9, 9, 9]
            else:
                # Start of a new pattern
                output[i] = [6, 9, 9]
                if i > 0 and output[i-1][1] != 4:
                    output[i-1] = [6, 6, 9]
    
    # Keep original values where no special patterns apply
    for i in range(height):
        if output[i][1] not in [4, 9]:
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