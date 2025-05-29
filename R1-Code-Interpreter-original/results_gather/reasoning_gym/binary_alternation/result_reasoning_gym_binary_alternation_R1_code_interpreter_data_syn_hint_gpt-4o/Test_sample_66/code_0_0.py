def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps for pattern starting with '0' (e.g., "010101...")
    swaps_0 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            swaps_0 += 1
    
    # Calculate swaps for pattern starting with '1' (e.g., "101010...")
    swaps_1 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            swaps_1 += 1
    
    # Each swap fixes two positions, so divide by 2
    if n % 2 == 0:
        # If length is even, both patterns are possible
        return min(swaps_0 // 2, swaps_1 // 2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_0 // 2
        else:
            return swaps_1 // 2

# Given binary string
binary_string = "001100010111110011010"
result = min_swaps_to_alternate(binary_string)
print(result)