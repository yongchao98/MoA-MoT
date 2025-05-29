def transform_grid(input_grid):
    # Convert input to list of integers
    nums = [int(x) for x in input_grid.split()]
    
    # Find the non-zero sequence
    non_zero_val = next(x for x in nums if x != 0)
    start_pos = next(i for i, x in enumerate(nums) if x != 0)
    sequence_length = len([x for x in nums if x != 0])
    
    # Create output grid with 4-position left shift
    new_start = max(0, start_pos - 4)  # Shift left by 4, but not below 0
    output = [0] * 29  # Initialize with zeros
    
    # Place the non-zero sequence
    for i in range(sequence_length):
        if new_start + i < 29:  # Ensure we don't exceed grid length
            output[new_start + i] = non_zero_val
    
    return ' '.join(map(str, output))

test_input = "0 0 0 0 0 7 7 7 7 7 7 7 7 7 7 7 7 7 0 0 0 0 0 0 0 0 0 0 0"
result = transform_grid(test_input)
print(result)