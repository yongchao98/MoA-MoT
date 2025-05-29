def analyze_and_shift(input_grid):
    # Convert string to list of integers
    grid = [int(x) for x in input_grid.split()]
    
    # Find the non-zero sequence
    non_zero = []
    start_idx = -1
    for i, val in enumerate(grid):
        if val != 0:
            if start_idx == -1:
                start_idx = i
            non_zero.append(val)
    
    # Create output grid
    output = [0] * len(grid)
    
    # Shift the non-zero sequence 4 positions to the right
    new_start = start_idx + 4
    for i, val in enumerate(non_zero):
        if new_start + i < len(grid):
            output[new_start + i] = val
    
    return ' '.join(map(str, output))

# Test input
test = "0 0 0 0 0 0 0 0 0 0 0 0 0 8 8 8 8 0 0 0 0 0 0 0 0 0 0"
print(analyze_and_shift(test))