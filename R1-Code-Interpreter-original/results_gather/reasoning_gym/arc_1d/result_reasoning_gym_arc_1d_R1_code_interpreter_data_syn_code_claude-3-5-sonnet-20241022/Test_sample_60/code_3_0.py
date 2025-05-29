def transform_grid(input_grid):
    output = [0] * len(input_grid)
    n = len(input_grid)
    
    # Find the end positions of consecutive sequences
    i = 0
    while i < n:
        if input_grid[i] != 0:
            # Find the end of this sequence
            start = i
            while i < n and input_grid[i] != 0:
                i += 1
            # i-1 is now the end of the sequence
            # Create pattern centered on the end of sequence
            pattern_center = i - 1
            pattern_start = max(0, pattern_center - 2)
            pattern_end = min(n, pattern_center + 3)
            
            # Fill the pattern
            for j in range(pattern_start, pattern_end):
                output[j] = input_grid[start]
        else:
            i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Print result
print(transform_grid(test_input))