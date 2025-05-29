def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to rearrange into an alternating pattern
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps for pattern "010101..."
    swaps_pattern_a = 0
    swaps_pattern_b = 0
    
    for i in range(n):
        expected_char_a = '0' if i % 2 == 0 else '1'
        expected_char_b = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_a:
            swaps_pattern_a += 1
        if s[i] != expected_char_b:
            swaps_pattern_b += 1
    
    # Since each swap fixes two positions, we divide by 2
    if n % 2 == 0:
        # If length is even, both patterns are possible
        return min(swaps_pattern_a // 2, swaps_pattern_b // 2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_pattern_a // 2
        else:
            return swaps_pattern_b // 2

# Given binary string
binary_string = "100111101001011000100"
result = min_swaps_to_alternate(binary_string)
print(result)