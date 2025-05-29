def transform_sequence(input_str):
    # Convert input to list of integers
    nums = [int(x) for x in input_str.split()]
    
    # Verify length
    if len(nums) != 29:
        raise ValueError("Input must have exactly 29 numbers")
    
    # Find the middle sequence
    main_num = nums[0]
    start_idx = None
    end_idx = None
    
    # Find start and end of middle sequence
    for i in range(len(nums)):
        if nums[i] != main_num and start_idx is None:
            start_idx = i
        elif nums[i] == main_num and start_idx is not None and end_idx is None:
            end_idx = i
            break
    
    if end_idx is None:
        end_idx = len(nums)
    
    # Create result array
    result = [main_num] * 29
    
    # Copy the middle sequence shifted right by 3
    middle_sequence = nums[start_idx:end_idx]
    new_start = start_idx + 3
    
    for i in range(len(middle_sequence)):
        result[new_start + i] = middle_sequence[i]
    
    return ' '.join(map(str, result))

# Test input
test_input = "1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1"
print("Final output:", transform_sequence(test_input))