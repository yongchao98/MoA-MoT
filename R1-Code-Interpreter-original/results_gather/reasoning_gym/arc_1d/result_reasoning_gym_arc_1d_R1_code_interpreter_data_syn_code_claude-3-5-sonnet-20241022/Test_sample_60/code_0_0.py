def transform_grid(input_grid):
    output = [0] * len(input_grid)
    
    # For each position in input
    for i in range(len(input_grid)):
        if input_grid[i] != 0:
            # Extend pattern 2 positions left and right
            for j in range(max(0, i-2), min(len(input_grid), i+3)):
                output[j] = input_grid[i]
    
    # Convert to string format
    return ' '.join(map(str, output))

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Print result
print(transform_grid(test_input))