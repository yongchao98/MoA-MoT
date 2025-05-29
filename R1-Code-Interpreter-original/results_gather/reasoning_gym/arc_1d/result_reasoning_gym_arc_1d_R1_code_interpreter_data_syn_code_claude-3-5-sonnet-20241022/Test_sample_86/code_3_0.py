def apply_rule(input_str):
    # Convert input to list of integers
    grid = [int(x) for x in input_str.split()]
    n = len(grid)
    result = [0] * n
    
    # Find non-zero sequences
    nonzero_nums = []
    for x in grid:
        if x != 0:
            nonzero_nums.append(x)
    
    # Find first non-zero position
    first_pos = next(i for i, x in enumerate(grid) if x != 0)
    
    # Rule: Move all numbers 4 positions to the left
    new_start = (first_pos - 4) % n
    
    # Place numbers
    for i, num in enumerate(nonzero_nums):
        result[(new_start + i) % n] = num
    
    return ' '.join(map(str, result))

# Test input
test_input = "8 1 2 0 0 0 0 0 0 0 8 7 2 8 4 5 1 8 2"
print(apply_rule(test_input))