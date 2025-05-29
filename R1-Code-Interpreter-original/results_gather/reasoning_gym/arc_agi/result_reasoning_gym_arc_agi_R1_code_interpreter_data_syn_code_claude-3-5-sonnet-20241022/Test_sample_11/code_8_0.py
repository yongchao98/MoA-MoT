def process_grid(grid):
    # Make a copy of the input grid using list comprehension
    result = [[x for x in row] for row in grid]
    
    # Get the dimensions
    rows, cols = len(grid), len(grid[0])
    
    # Process each position
    for i in range(rows):
        for j in range(cols):
            # Skip if not zero
            if grid[i][j] != 0:
                continue
                
            # Get the replacement value based on position
            if i <= 4:  # Upper half (including middle)
                replacement = grid[0][j]  # First row value
            else:  # Lower half
                replacement = grid[-1][j]  # Last row value
                
            # Replace the zero
            result[i][j] = replacement
    
    # Print for verification
    for row in result:
        print(*row, sep=' ')
    
    return result

# Test case
input_grid = [
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,0,8,8,8,0,8,8,6],
    [4,8,8,8,0,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,0,8,8,6],
    [4,8,8,0,8,8,8,8,8,6],
    [4,8,8,8,8,8,8,8,8,6],
    [4,8,8,8,8,0,8,8,0,6],
    [4,0,8,8,8,8,8,8,8,6]
]

process_grid(input_grid)