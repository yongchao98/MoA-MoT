def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
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
    
    # If count_0 > count_1, we can only use pattern "010101..."
    if count_0 > count_1:
        return swaps_010 // 2
    # If count_1 > count_0, we can only use pattern "101010..."
    elif count_1 > count_0:
        return swaps_101 // 2
    # If count_0 == count_1, we can use either pattern
    else:
        return min(swaps_010, swaps_101) // 2

# Given binary string
binary_string = "100101010010011010111001010"
result = min_swaps_to_alternate(binary_string)
print(result)