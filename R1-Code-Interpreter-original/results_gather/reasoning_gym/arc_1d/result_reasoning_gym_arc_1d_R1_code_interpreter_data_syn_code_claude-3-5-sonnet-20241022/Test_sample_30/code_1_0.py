def apply_rule(input_str):
    nums = [int(x) for x in input_str.split()]
    
    # Find the non-zero sequence
    non_zeros = [(i, num) for i, num in enumerate(nums) if num != 0]
    sequence = [num for _, num in non_zeros]
    
    result = [0] * 12  # Initialize output with zeros
    
    if len(sequence) == 4:
        # Split into two pairs
        result[0] = sequence[2]  # First pair goes to start
        result[1] = sequence[3]
        result[10] = sequence[0]  # Second pair goes to end
        result[11] = sequence[1]
    elif len(sequence) == 3:
        # Place in middle starting at position 5
        result[5] = sequence[0]
        result[6] = sequence[1]
        result[7] = sequence[2]
    elif len(sequence) == 2:
        # Place in middle starting at position 3
        result[3] = sequence[0]
        result[4] = sequence[1]
    
    return result

# Test with the new input
test_input = "6 3 4 2 6 9 8 1 7 0 0 0"
result = apply_rule(test_input)
print(" ".join(map(str, result)))