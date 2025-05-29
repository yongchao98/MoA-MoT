def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern starting with '0' (i.e., "010101...")
    swaps_0 = 0
    swaps_1 = 0
    for i in range(n):
        expected_char_0 = '0' if i % 2 == 0 else '1'
        expected_char_1 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_0:
            swaps_0 += 1
        if s[i] != expected_char_1:
            swaps_1 += 1
    
    # Since each swap fixes two positions, we divide by 2
    if n % 2 == 0:
        return min(swaps_0, swaps_1) // 2
    else:
        # If the string length is odd, we need to consider the valid pattern
        if count_0 > count_1:
            return swaps_0 // 2
        else:
            return swaps_1 // 2

# Given binary string
binary_string = "1111111000010000"
result = min_swaps_to_alternate(binary_string)
print(result)