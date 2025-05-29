def solve_grid(input_str):
    # Convert input to list
    grid = [int(x) for x in input_str.split()]
    
    # Find the majority number
    from collections import Counter
    majority = Counter(grid).most_common(1)[0][0]
    
    # Find the first and last position of the non-majority sequence
    start = None
    end = None
    for i, num in enumerate(grid):
        if num != majority:
            if start is None:
                start = i
            end = i
    
    # If no different sequence found, return input as is
    if start is None:
        return input_str
    
    # Create output grid
    output = [majority] * len(grid)
    
    # Move sequence 2 positions left
    sequence = grid[start:end+1]
    new_start = max(0, start - 2)
    for i, num in enumerate(sequence):
        if new_start + i < len(grid):
            output[new_start + i] = num
    
    return ' '.join(map(str, output))

# Test with the given input
test_input = "0 0 0 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0"
print(solve_grid(test_input))