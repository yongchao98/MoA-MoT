def transform_final(input_str):
    nums = [int(x) for x in input_str.split()]
    main_num = nums[0]
    result = [main_num] * 29
    
    # Find the middle sequence
    start_idx = None
    end_idx = None
    for i in range(len(nums)):
        if nums[i] != main_num and start_idx is None:
            start_idx = i
        elif nums[i] == main_num and start_idx is not None and end_idx is None:
            end_idx = i
            break
    
    if end_idx is None:
        end_idx = len(nums)
    
    # The middle sequence shifts right by 3
    middle_sequence = nums[start_idx:end_idx]
    new_start = start_idx + 3
    
    # Place the middle sequence
    for i in range(len(middle_sequence)):
        result[new_start + i] = middle_sequence[i]
    
    return ' '.join(map(str, result))

test_input = "1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1"
print("Final correct output:", transform_final(test_input))