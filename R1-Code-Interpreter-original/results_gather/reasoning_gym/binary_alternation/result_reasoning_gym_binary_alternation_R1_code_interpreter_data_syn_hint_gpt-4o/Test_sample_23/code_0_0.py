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
    
    # Calculate minimum swaps needed
    if n % 2 == 0:
        # If length is even, both patterns are possible
        min_swaps = min(mismatches_pattern1 // 2, mismatches_pattern2 // 2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            min_swaps = mismatches_pattern1 // 2
        else:
            min_swaps = mismatches_pattern2 // 2
    
    return min_swaps

# Given binary string
binary_string = "10001111001110"
print(min_swaps_to_alternate(binary_string))