def min_swaps_to_alternate(s):
    n = len(s)
    pattern1_mismatches = 0
    pattern2_mismatches = 0
    
    # Calculate mismatches for both patterns
    for i in range(n):
        expected_char_pattern1 = '0' if i % 2 == 0 else '1'
        expected_char_pattern2 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_pattern1:
            pattern1_mismatches += 1
        if s[i] != expected_char_pattern2:
            pattern2_mismatches += 1
    
    # Check if it's possible to make the string alternating
    if pattern1_mismatches % 2 != 0 or pattern2_mismatches % 2 != 0:
        return -1
    
    # Calculate the minimum swaps needed
    min_swaps = min(pattern1_mismatches, pattern2_mismatches) // 2
    return min_swaps

# Given binary string
binary_string = "01000111101110"
result = min_swaps_to_alternate(binary_string)
print(result)