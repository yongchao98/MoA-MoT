def transform_grid(input_grid):
    output = [0] * len(input_grid)
    n = len(input_grid)
    
    i = 0
    last_pattern_end = 0  # Keep track of where the last pattern ended
    
    while i < n:
        if input_grid[i] != 0:
            # Find the sequence length
            seq_start = i
            while i < n and input_grid[i] != 0:
                i += 1
            seq_length = i - seq_start
            
            # Calculate middle position of sequence
            middle = seq_start + (seq_length - 1) // 2
            
            # Create pattern centered on middle position
            pattern_start = max(0, middle - 2)
            pattern_end = min(n, pattern_start + 5)  # Ensure exactly 5 positions
            
            # Only apply pattern if it doesn't overlap with previous pattern
            if pattern_start >= last_pattern_end:
                for j in range(pattern_start, pattern_end):
                    output[j] = input_grid[seq_start]
                last_pattern_end = pattern_end
        else:
            i += 1
    
    return ' '.join(map(str, output))

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Print result
print(transform_grid(test_input))