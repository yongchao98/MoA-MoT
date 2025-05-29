def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to rearrange into an alternating pattern
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
    if n % 2 == 0:
        # If length is even, both patterns are possible
        min_swaps = min(mismatches_01 // 2, mismatches_10 // 2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            min_swaps = mismatches_01 // 2
        else:
            min_swaps = mismatches_10 // 2
    
    return min_swaps

# Given binary string
binary_string = "010101000000011111111100001"
print(min_swaps_to_alternate(binary_string))