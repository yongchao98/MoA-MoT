def transform_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Initialize output with first three columns
    output = [[rows[i][j] for j in range(3)] for i in range(height)]
    
    # Process patterns from top to bottom
    for i in range(height):
        # Check for 4 pattern
        if output[i][1] == 4:
            # Affect current row and next two rows only
            for j in range(i, min(i+3, height)):
                output[j] = [4, 4, 4]
            # Skip the affected rows
            i += 2
            continue
            
        # Check for 9 pattern (only if not in a 4 pattern)
        if output[i][1] == 9 and output[i][0] != 4:
            # Set current row to [9,9,9]
            output[i] = [9, 9, 9]
            # Set row above to [6,9,9] if exists and not in 4 pattern
            if i > 0 and output[i-1][1] != 4:
                output[i-1] = [6, 9, 9]
            # Set two rows above to [6,6,9] if exists and not in 4 pattern
            if i > 1 and output[i-2][1] != 4:
                output[i-2] = [6, 6, 9]
    
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