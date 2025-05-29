def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern starting with '0' (e.g., "010101...")
    swaps_start_with_0 = 0
    swaps_start_with_1 = 0
    
    for i in range(n):
        expected_char_0 = '0' if i % 2 == 0 else '1'
        expected_char_1 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_0:
            swaps_start_with_0 += 1
        if s[i] != expected_char_1:
            swaps_start_with_1 += 1
    
    # Calculate the minimum swaps needed
    if count_0 == count_1:
        # Both patterns are possible
        return min(swaps_start_with_0 // 2, swaps_start_with_1 // 2)
    elif count_0 > count_1:
        # Only pattern starting with '0' is possible
        return swaps_start_with_0 // 2
    else:
        # Only pattern starting with '1' is possible
        return swaps_start_with_1 // 2

# Given binary string
binary_string = "00111110010111001000001"
result = min_swaps_to_alternate(binary_string)
print(result)