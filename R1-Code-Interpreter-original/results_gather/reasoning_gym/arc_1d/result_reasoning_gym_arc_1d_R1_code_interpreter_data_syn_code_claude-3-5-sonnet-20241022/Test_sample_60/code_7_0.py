def transform_grid(input_grid):
    output = [0] * len(input_grid)
    n = len(input_grid)
    
    i = 0
    while i < n:
        if input_grid[i] != 0:
            # Start of a sequence - create a pattern right here
            pattern_start = i
            pattern_end = min(n, pattern_start + 5)  # Exactly 5 positions
            
            # Fill the pattern
            for j in range(pattern_start, pattern_end):
                output[j] = input_grid[i]
            
            # Skip to end of sequence
            while i < n and input_grid[i] != 0:
                i += 1
        else:
            i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Print result
print(transform_grid(test_input))