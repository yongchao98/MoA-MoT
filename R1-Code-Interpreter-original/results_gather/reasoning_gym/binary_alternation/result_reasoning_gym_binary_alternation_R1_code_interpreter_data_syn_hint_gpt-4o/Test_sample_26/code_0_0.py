def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check if it's possible to make the string alternating
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern "010101..."
    swaps_010 = 0
    swaps_101 = 0
    
    for i in range(n):
        expected_char_010 = '0' if i % 2 == 0 else '1'
        expected_char_101 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_010:
            swaps_010 += 1
        if s[i] != expected_char_101:
            swaps_101 += 1
    
    # Since each swap fixes two positions, divide by 2
    return min(swaps_010, swaps_101) // 2

# Given binary string
binary_string = "0101010011110000101111001100"
result = min_swaps_to_alternate(binary_string)
print(result)