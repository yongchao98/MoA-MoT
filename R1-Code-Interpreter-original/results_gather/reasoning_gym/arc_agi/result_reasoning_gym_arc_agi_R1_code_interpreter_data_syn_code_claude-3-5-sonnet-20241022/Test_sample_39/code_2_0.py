def transform_grid(input_grid):
    # Convert input strings to 2D array
    grid = [row.split() for row in input_grid]
    output = [row[:] for row in grid]
    
    # Pattern repeats every 4 rows
    for i in range(0, len(grid), 4):
        group = grid[i:i+4]
        
        # First row in group
        if i < len(grid):
            if grid[i].count('0') == 1:
                pos = grid[i].index('0')
                output[i][pos] = '5'
                
        # Second row in group
        if i+1 < len(grid):
            continue  # Keep original
            
        # Third row in group
        if i+2 < len(grid):
            if '0' in grid[i+2]:
                pos = grid[i+2].index('0')
                if pos == 2:  # If '0' is in last position
                    output[i+2][pos] = '5'
                
        # Fourth row in group
        if i+3 < len(grid):
            continue  # Keep original
    
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