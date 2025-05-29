def process_grid(input_grid):
    # Convert input to 2D array
    rows = [[int(x) for x in row.split()] for row in input_grid.strip().split('\n')]
    height = len(rows)
    
    # Extract groups (excluding '5's)
    groups = []
    for col in range(0, len(rows[0]), 4):
        if col + 3 <= len(rows[0]):
            group = [[rows[i][col], rows[i][col+1], rows[i][col+2]] for i in range(height)]
            groups.append(group)
    
    # Process propagation from right to left
    output = [[6, 6, 6] for _ in range(height)]  # Default state
    
    # Process each group from right to left
    for group in reversed(groups):
        for row in range(height):
            for col in range(3):
                if group[row][col] in [9, 7, 4]:  # Special numbers that can propagate
                    # Update output state
                    output[row][col] = group[row][col]
    
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

result = process_grid(test_input)
print("Predicted output:")
print(result)