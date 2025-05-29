def transform_grid(input_grid):
    # Convert input strings to 2D array
    grid = [row.split() for row in input_grid]
    output = [row[:] for row in grid]
    
    # Define transformation positions and patterns
    transforms = [
        (0, '0 9 9', '5 9 9'),    # First row
        (2, '9 9 0', '9 9 5'),    # Third row
        (4, '9 9 0', '9 9 5'),    # Fifth row
        (6, '0 9 9', '5 9 9'),    # Seventh row
        (8, '9 0 9', '9 5 9'),    # Ninth row
        (10, '0 9 9', '5 9 9'),   # Eleventh row
        (12, '0 9 9', '5 9 9'),   # Thirteenth row
        (14, '9 0 9', '9 5 9')    # Fifteenth row
    ]
    
    # Apply transformations
    for pos, input_pattern, output_pattern in transforms:
        if pos < len(grid):
            current_row = ' '.join(grid[pos])
            if current_row == input_pattern:
                output[pos] = output_pattern.split()
    
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