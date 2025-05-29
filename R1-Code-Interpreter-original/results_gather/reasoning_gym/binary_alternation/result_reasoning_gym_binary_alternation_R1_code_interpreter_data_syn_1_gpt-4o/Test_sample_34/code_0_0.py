def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for pattern starting with '0' (e.g., "010101...")
    mismatches_0 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            mismatches_0 += 1
    
    # Calculate mismatches for pattern starting with '1' (e.g., "101010...")
    mismatches_1 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            mismatches_1 += 1
    
    # Calculate swaps needed
    swaps_0 = mismatches_0 // 2
    swaps_1 = mismatches_1 // 2
    
    # If the string length is even, both patterns are possible
    if n % 2 == 0:
        return min(swaps_0, swaps_1)
    else:
        # If the string length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_0
        else:
            return swaps_1

# Given binary string
binary_string = "011000110101"
result = min_swaps_to_alternate(binary_string)
print(result)