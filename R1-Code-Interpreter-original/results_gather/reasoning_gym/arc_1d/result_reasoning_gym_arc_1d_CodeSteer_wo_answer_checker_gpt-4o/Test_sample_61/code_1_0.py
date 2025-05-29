def transform_grid(input_grid):
    output_grid = input_grid[:]
    n = len(input_grid)
    
    for i in range(n):
        # Use a sliding window to examine the current element and its neighbors
        left = input_grid[i-1] if i > 0 else None
        current = input_grid[i]
        right = input_grid[i+1] if i < n-1 else None
        
        # Apply transformation rules based on the pattern
        if current == 2:
            if (left in [0, 2, None]) and (right in [2, None]):
                output_grid[i] = 8
            elif left == 8:
                output_grid[i] = 8
    
    return output_grid

# Test input
test_input = [2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
test_output = transform_grid(test_input)

# Print the output in the required format
print("<<<" + " ".join(map(str, test_output)) + ">>>")