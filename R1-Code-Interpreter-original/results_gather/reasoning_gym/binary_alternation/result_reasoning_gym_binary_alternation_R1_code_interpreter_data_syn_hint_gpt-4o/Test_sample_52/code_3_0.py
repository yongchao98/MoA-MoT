def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to form an alternating string
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for both patterns
    mismatch_pattern1 = 0  # "010101..."
    mismatch_pattern2 = 0  # "101010..."
    
    for i in range(n):
        expected_char_pattern1 = '0' if i % 2 == 0 else '1'
        expected_char_pattern2 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_pattern1:
            mismatch_pattern1 += 1
        if s[i] != expected_char_pattern2:
            mismatch_pattern2 += 1
    
    # Determine the minimum swaps needed
    if n % 2 == 0:
        # If length is even, both patterns are possible
        min_swaps = min(mismatch_pattern1 // 2, mismatch_pattern2 // 2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            min_swaps = mismatch_pattern1 // 2
        else:
            min_swaps = mismatch_pattern2 // 2
    
    return min_swaps

# Given binary string
binary_string = "10000101111010101101"
result = min_swaps_to_alternate(binary_string)
print(result)