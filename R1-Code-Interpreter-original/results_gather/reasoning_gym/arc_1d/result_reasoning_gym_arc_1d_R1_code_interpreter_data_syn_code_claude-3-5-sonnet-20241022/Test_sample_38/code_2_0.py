def transform_sequence(input_str):
    # Convert input to list
    nums = [int(x) for x in input_str.split()]
    
    # Count leading zeros
    leading_zeros = 0
    for n in nums:
        if n == 0:
            leading_zeros += 1
        else:
            break
    
    # Find last non-zero number
    last_non_zero = len(nums) - 1
    while last_non_zero >= 0 and nums[last_non_zero] == 0:
        last_non_zero -= 1
    
    # Create output list
    output = []
    
    # Add leading zeros if more than 4
    if leading_zeros > 4:
        output.extend([0] * leading_zeros)
    
    # Add the sequence from first non-zero to last non-zero
    for i in range(leading_zeros, last_non_zero + 1):
        output.append(nums[i])
    
    # Fill remaining positions with zeros
    while len(output) < len(nums):
        output.append(0)
    
    print(' '.join(map(str, output)))

# Test input
test = "0 0 0 0 6 3 5 1 7 9 1 2 4 9 6 3 2 0 0"
transform_sequence(test)