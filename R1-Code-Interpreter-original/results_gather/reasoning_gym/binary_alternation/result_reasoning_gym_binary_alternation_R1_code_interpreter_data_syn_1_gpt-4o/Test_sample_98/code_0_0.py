def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Generate the two possible alternating patterns
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Count mismatches for both patterns
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # Calculate the minimum swaps needed
    if count_0 == count_1:
        # If counts are equal, we can use either pattern
        return min(mismatches1, mismatches2) // 2
    elif count_0 > count_1:
        # If more 0s, pattern1 must be used
        return mismatches1 // 2
    else:
        # If more 1s, pattern2 must be used
        return mismatches2 // 2

# Given binary string
binary_string = "00111100111010101111000000"
result = min_swaps_to_alternate(binary_string)
print(result)