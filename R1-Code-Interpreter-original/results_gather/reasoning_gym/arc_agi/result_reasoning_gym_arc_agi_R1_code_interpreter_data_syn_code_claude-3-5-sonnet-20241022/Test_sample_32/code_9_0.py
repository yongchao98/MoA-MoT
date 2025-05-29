def transform_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Initialize output with first three columns
    output = [[rows[i][j] for j in range(3)] for i in range(height)]
    
    # Process from top to bottom
    i = 0
    while i < height:
        # Handle 4 pattern
        if output[i][1] == 4:
            # Current row keeps its values
            # Next two rows become [4,4,4]
            if i + 1 < height:
                output[i+1] = [4, 4, 4]
            if i + 2 < height:
                output[i+2] = [4, 4, 4]
            i += 3
            continue
        
        # Handle 9 pattern
        if output[i][1] == 9:
            # Current row becomes [9,9,9]
            output[i] = [9, 9, 9]
            # Previous row becomes [6,9,9] if it exists and isn't part of a 4 pattern
            if i > 0 and output[i-1][1] != 4:
                output[i-1] = [6, 9, 9]
            # Two rows up becomes [6,6,9] if it exists and isn't part of a 4 pattern
            if i > 1 and output[i-2][1] != 4:
                output[i-2] = [6, 6, 9]
        i += 1
    
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