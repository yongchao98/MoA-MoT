def transform_grid(input_grid):
    # Convert input strings to 2D array
    grid = [row.split() for row in input_grid]
    output = [row[:] for row in grid]
    
    # Process each row
    for i in range(len(grid)):
        row = grid[i]
        # Pattern 1: "0 9 9" -> "5 9 9"
        if row == ['0', '9', '9'] and (i % 4 == 0 or i % 6 == 2):
            output[i] = ['5', '9', '9']
        # Pattern 2: "9 9 0" -> "9 9 5"
        elif row == ['9', '9', '0'] and (i % 4 == 2):
            output[i] = ['9', '9', '5']
        # Pattern 3: "9 0 9" -> "9 5 9"
        elif row == ['9', '0', '9'] and (i % 8 == 4):
            output[i] = ['9', '5', '9']
    
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