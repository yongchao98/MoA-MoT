def find_middle_sequence(nums):
    # Find the different number in the middle
    main_num = nums[0]
    start_idx = None
    end_idx = None
    
    # Find start of different sequence
    for i in range(len(nums)):
        if nums[i] != main_num:
            start_idx = i
            break
    
    # Find end of different sequence
    for i in range(start_idx, len(nums)):
        if nums[i] == main_num:
            end_idx = i
            break
    if end_idx is None:
        end_idx = len(nums)
    
    return start_idx, end_idx

def transform_grid(input_str):
    # Convert input string to list of integers
    nums = [int(x) for x in input_str.split()]
    result = nums.copy()
    
    # Find the middle sequence
    start_idx, end_idx = find_middle_sequence(nums)
    seq_length = end_idx - start_idx
    
    # Shift the middle sequence right by 3 positions
    main_num = nums[0]
    diff_num = nums[start_idx]
    
    # Reset the original positions to main number
    for i in range(start_idx, end_idx):
        result[i] = main_num
    
    # Place the sequence in new position (shifted right by 3)
    for i in range(seq_length):
        result[start_idx + 3 + i] = diff_num
    
    return ' '.join(map(str, result))

# Test input
test_input = "1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1"
print("Transformed output:", transform_grid(test_input))