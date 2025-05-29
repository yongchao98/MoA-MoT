def apply_correct_rule(input_grid):
    # Convert input string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Find position of '2'
    marker_pos = grid.index(2)
    
    # Find the number that follows '2'
    if marker_pos < len(grid) - 1:
        replacement = grid[marker_pos + 1]
    
    # Keep everything before marker_pos
    # Replace everything from marker_pos onwards with the replacement number
    output = grid[:marker_pos] + [replacement] * (len(grid) - marker_pos)
    
    # Convert back to string
    return ' '.join(map(str, output))

# Test input
test = "2 0 0 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4"
result = apply_correct_rule(test)
print(result)