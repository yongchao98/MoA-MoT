def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between counts of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for both patterns
    mismatches_01 = 0  # Pattern starting with '0'
    mismatches_10 = 0  # Pattern starting with '1'
    
    for i in range(n):
        expected_char_01 = '0' if i % 2 == 0 else '1'
        expected_char_10 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_01:
            mismatches_01 += 1
        if s[i] != expected_char_10:
            mismatches_10 += 1
    
    # Calculate swaps needed
    swaps_01 = mismatches_01 // 2
    swaps_10 = mismatches_10 // 2
    
    # Return the minimum swaps needed
    return min(swaps_01, swaps_10)

# Given binary string
binary_string = "1100101011000010110101011"
result = min_swaps_to_alternate(binary_string)
print(result)