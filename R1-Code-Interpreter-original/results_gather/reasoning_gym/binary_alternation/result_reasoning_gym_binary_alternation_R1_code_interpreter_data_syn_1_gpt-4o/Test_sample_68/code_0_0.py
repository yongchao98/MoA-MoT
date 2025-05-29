def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for both patterns
    mismatches_pattern1 = 0  # Pattern "010101..."
    mismatches_pattern2 = 0  # Pattern "101010..."
    
    for i in range(n):
        expected_char_pattern1 = '0' if i % 2 == 0 else '1'
        expected_char_pattern2 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_pattern1:
            mismatches_pattern1 += 1
        if s[i] != expected_char_pattern2:
            mismatches_pattern2 += 1
    
    # Calculate swaps needed
    swaps_pattern1 = mismatches_pattern1 // 2
    swaps_pattern2 = mismatches_pattern2 // 2
    
    # Return the minimum swaps needed
    if count_0 == count_1:
        return min(swaps_pattern1, swaps_pattern2)
    elif count_0 > count_1:
        return swaps_pattern1
    else:
        return swaps_pattern2

# Given binary string
binary_string = "00010111011010"
result = min_swaps_to_alternate(binary_string)
print(result)