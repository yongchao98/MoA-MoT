def transform_grid(input_grid):
    # Convert input strings to 2D array
    grid = [row.split() for row in input_grid]
    output = [row[:] for row in grid]
    
    # Define the pattern of transformations
    # Each tuple contains (row_offset, position_to_check, position_to_transform)
    patterns = [
        (0, 1, 1),    # middle position in first row
        (2, 2, 2),    # right position in third row
        (4, 2, 2),    # right position in fifth row
        (6, 0, 0),    # left position in seventh row
        (8, 1, 1),    # middle position in ninth row
        (10, 0, 0),   # left position in eleventh row
        (12, 0, 0)    # left position in thirteenth row
    ]
    
    # Apply transformations
    for base_row in range(0, len(grid), 14):  # Pattern repeats every 14 rows
        for row_offset, check_pos, transform_pos in patterns:
            current_row = base_row + row_offset
            if current_row < len(grid):
                if grid[current_row][check_pos] == '0':
                    output[current_row][transform_pos] = '5'
    
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