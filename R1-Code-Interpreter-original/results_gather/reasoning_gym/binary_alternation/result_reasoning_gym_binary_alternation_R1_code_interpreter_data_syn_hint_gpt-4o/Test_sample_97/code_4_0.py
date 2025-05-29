def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Pattern 1: Starts with '0'
    pattern1_mismatches = 0
    # Pattern 2: Starts with '1'
    pattern2_mismatches = 0
    
    for i in range(n):
        expected_char_pattern1 = '0' if i % 2 == 0 else '1'
        expected_char_pattern2 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_pattern1:
            pattern1_mismatches += 1
        if s[i] != expected_char_pattern2:
            pattern2_mismatches += 1
    
    # Calculate swaps needed
    swaps_pattern1 = pattern1_mismatches // 2
    swaps_pattern2 = pattern2_mismatches // 2
    
    # If counts of 0s and 1s are equal, both patterns are possible
    if count_0 == count_1:
        return min(swaps_pattern1, swaps_pattern2)
    # If there are more 0s, pattern starting with '0' is the only possible
    elif count_0 > count_1:
        return swaps_pattern1
    # If there are more 1s, pattern starting with '1' is the only possible
    else:
        return swaps_pattern2

# Given binary string
binary_string = "110010100010110001010101101001"
result = min_swaps_to_alternate(binary_string)
print(result)