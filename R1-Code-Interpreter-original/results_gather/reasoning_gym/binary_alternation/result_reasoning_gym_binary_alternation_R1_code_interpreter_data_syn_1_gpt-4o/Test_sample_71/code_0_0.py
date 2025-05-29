def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: "010101..."
    mismatches_01 = 0
    mismatches_10 = 0
    
    for i in range(n):
        expected_char_01 = '0' if i % 2 == 0 else '1'
        expected_char_10 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_01:
            mismatches_01 += 1
        if s[i] != expected_char_10:
            mismatches_10 += 1
    
    # Check if conversion is possible
    if mismatches_01 % 2 != 0 and mismatches_10 % 2 != 0:
        return -1
    
    # Calculate minimum swaps
    swaps_01 = mismatches_01 // 2
    swaps_10 = mismatches_10 // 2
    
    if mismatches_01 % 2 == 0 and mismatches_10 % 2 == 0:
        return min(swaps_01, swaps_10)
    elif mismatches_01 % 2 == 0:
        return swaps_01
    else:
        return swaps_10

# Given binary string
binary_string = "10101110001100100100010101"
result = min_swaps_to_alternate(binary_string)
print(result)