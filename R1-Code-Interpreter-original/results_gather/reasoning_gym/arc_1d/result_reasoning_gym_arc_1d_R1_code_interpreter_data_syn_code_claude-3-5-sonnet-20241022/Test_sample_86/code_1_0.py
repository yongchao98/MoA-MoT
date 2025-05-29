def apply_rule(input_str):
    # Convert input to list of integers
    grid = [int(x) for x in input_str.split()]
    n = len(grid)
    result = [0] * n
    
    # Find all non-zero numbers in sequence
    nonzero_nums = [x for x in grid if x != 0]
    
    # If there are non-zero numbers, place them starting at the appropriate position
    if nonzero_nums:
        # Find the position of the first non-zero number in input
        first_nonzero_pos = next(i for i, x in enumerate(grid) if x != 0)
        # Calculate new starting position (4 positions to the right)
        new_start = (first_nonzero_pos + 4) % n
        
        # Place the numbers
        for i, num in enumerate(nonzero_nums):
            result[(new_start + i) % n] = num
    
    return ' '.join(map(str, result))

# Test input
test_input = "8 1 2 0 0 0 0 0 0 0 8 7 2 8 4 5 1 8 2"
print(apply_rule(test_input))