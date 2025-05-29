def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for both patterns
    mismatches_01 = 0  # Pattern starting with '0'
    mismatches_10 = 0  # Pattern starting with '1'
    
    for i in range(n):
        expected_char_01 = '0' if i % 2 == 0 else '1'
        expected_char_10 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_01:
            mismatches_01 += 1
        if s[i] != expected_char_10:
            mismatches_10 += 1
    
    # Calculate minimum swaps
    if count_0 == count_1:
        # If counts are equal, we can choose either pattern
        min_swaps = min(mismatches_01, mismatches_10) // 2
    elif count_0 > count_1:
        # If more '0's, pattern must start with '0'
        min_swaps = mismatches_01 // 2
    else:
        # If more '1's, pattern must start with '1'
        min_swaps = mismatches_10 // 2
    
    return min_swaps

# Given binary string
binary_string = "010001101010011111100010100"
result = min_swaps_to_alternate(binary_string)
print(result)