def transform_grid(input_grid):
    output = [0] * len(input_grid)
    patterns = []
    
    # First collect all patterns
    for i in range(len(input_grid)):
        if input_grid[i] != 0:
            # Calculate the start and end positions for the 5-number pattern
            start = max(0, i - 2)
            end = min(len(input_grid), start + 5)
            # If we're near the right edge, adjust start position
            if end - start < 5:
                start = max(0, end - 5)
            # Store pattern info
            patterns.append((start, end, input_grid[i]))
    
    # Sort patterns by start position
    patterns.sort()
    
    # Apply patterns sequentially, ensuring no overlap
    current_end = 0
    for start, end, value in patterns:
        if start >= current_end:  # Only apply if not overlapping with previous pattern
            for j in range(start, end):
                output[j] = value
            current_end = end
        elif start < current_end and end > current_end:  # Partial overlap
            for j in range(current_end, end):
                output[j] = value
            current_end = end
    
    return ' '.join(map(str, output))

# Test input
test_input = [3, 3, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]

# Print result
print(transform_grid(test_input))