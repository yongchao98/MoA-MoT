def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between the number of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
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
    
    # Calculate swaps needed for each pattern
    swaps_01 = mismatches_01 // 2
    swaps_10 = mismatches_10 // 2
    
    # Return the minimum swaps needed
    if count_0 == count_1:
        return min(swaps_01, swaps_10)
    elif count_0 > count_1:
        return swaps_01
    else:
        return swaps_10

# Given binary string
binary_string = "101110110000100001"
result = min_swaps_to_alternate(binary_string)
print(result)