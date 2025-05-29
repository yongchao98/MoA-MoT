def find_pattern(input_str):
    # Convert string to list
    nums = input_str.split()
    
    # Find the longest sequence of zeros and its position
    max_zero_seq = 0
    max_zero_start = -1
    current_seq = 0
    
    for i in range(len(nums)):
        if nums[i] == '0':
            current_seq += 1
            if current_seq > max_zero_seq:
                max_zero_seq = current_seq
                max_zero_start = i - current_seq + 1
        else:
            current_seq = 0
    
    # Create output by shifting the zero sequence right
    output = nums.copy()
    non_zero = [x for x in nums if x != '0']
    zero_count = nums.count('0')
    
    # Reconstruct the output
    result = []
    zero_placed = 0
    non_zero_idx = 0
    
    for i in range(len(nums)):
        if i < max_zero_start + 1 or i >= max_zero_start + zero_count + 1:
            if non_zero_idx < len(non_zero):
                result.append(non_zero[non_zero_idx])
                non_zero_idx += 1
        else:
            result.append('0')
            zero_placed += 1
    
    print(f"Input: {nums}")
    print(f"Output: {result}")
    return ' '.join(result)

# Test with the given input
test_input = "6 6 6 0 0 0 0 0 0 0 0 6"
output = find_pattern(test_input)