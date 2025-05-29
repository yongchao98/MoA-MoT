def transform_grid(input_grid):
    output = [0] * len(input_grid)
    n = len(input_grid)
    
    # Find the start positions of consecutive sequences
    i = 0
    while i < n:
        if input_grid[i] != 0:
            # This is the start of a sequence
            value = input_grid[i]
            # Create pattern centered on this position
            pattern_start = max(0, i - 2)
            pattern_end = min(n, i + 3)
            
            # Fill the pattern
            for j in range(pattern_start, pattern_end):
                output[j] = value
                
            # Skip the rest of this sequence
            while i < n and input_grid[i] != 0:
                i += 1
        else:
            i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Print result
print(transform_grid(test_input))