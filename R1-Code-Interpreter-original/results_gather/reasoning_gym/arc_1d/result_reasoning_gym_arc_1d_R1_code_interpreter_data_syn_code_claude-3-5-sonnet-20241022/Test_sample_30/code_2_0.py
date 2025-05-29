def find_rule(input_str):
    nums = [int(x) for x in input_str.split()]
    result = [0] * 12
    
    # Find the first non-zero sequence
    start = -1
    end = -1
    for i, num in enumerate(nums):
        if num != 0:
            if start == -1:
                start = i
            end = i
        elif start != -1 and num == 0:
            break
    
    if start == -1:
        return result
        
    sequence = nums[start:end+1]
    
    # Apply the rules
    if len(sequence) == 4:
        # If sequence has 4 numbers and is near the end (like Example 1)
        if start >= 7:  # sequence starts late in the array
            # First two numbers of output are last two numbers of input sequence
            result[0] = sequence[2]
            result[1] = sequence[3]
            # Last two numbers of output are first two numbers of input sequence
            result[10] = sequence[0]
            result[11] = sequence[1]
    else:
        # For other cases, move the sequence 3 positions to the right
        new_start = start + 3
        for i, num in enumerate(sequence):
            result[new_start + i] = num
            
    return result

# Test with the new input
test_input = "6 3 4 2 6 9 8 1 7 0 0 0"
result = find_rule(test_input)
print(" ".join(map(str, result)))