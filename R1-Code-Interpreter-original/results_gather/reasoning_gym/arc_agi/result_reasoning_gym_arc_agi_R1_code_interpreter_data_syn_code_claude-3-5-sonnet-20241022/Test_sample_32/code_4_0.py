def transform_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Initialize output with first three columns
    output = [[rows[i][j] for j in range(3)] for i in range(height)]
    
    # Process from bottom to top
    for i in range(height-1, -1, -1):
        # Check for 4 pattern
        if output[i][1] == 4:
            # Affect current row and two rows above
            for j in range(max(0, i-2), i+1):
                output[j] = [4, 4, 4]
        
        # Check for 9 pattern
        elif output[i][1] == 9:
            # Current row becomes [9,9,9]
            output[i] = [9, 9, 9]
            # Row above becomes [6,9,9] if not already affected
            if i > 0 and output[i-1][1] != 4:
                output[i-1] = [6, 9, 9]
            # Two rows above becomes [6,6,9] if not already affected
            if i > 1 and output[i-2][1] != 4:
                output[i-2] = [6, 6, 9]
    
    # Convert to string format
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