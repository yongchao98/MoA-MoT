def transform_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Get the first group (first three columns)
    base = [[rows[i][j] for j in range(3)] for i in range(height)]
    
    # Process special patterns
    output = [row[:] for row in base]  # Start with base pattern
    
    # Process patterns from bottom to top
    for i in range(height-1, -1, -1):
        # Check for special number patterns
        if i < height-1:  # Look at relationships with row below
            if base[i][1] == 4:  # Pattern with 4
                for j in range(i-2, i+1):  # Affect 3 rows up
                    if 0 <= j < height:
                        output[j] = [4, 4, 4]
            elif base[i][1] == 9:  # Pattern with 9
                output[i] = [6, 9, 9]
                if i > 0:
                    output[i-1] = [6, 6, 9]
    
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