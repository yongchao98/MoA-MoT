def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Generate two possible alternating patterns
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Count mismatches for both patterns
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # Calculate minimum swaps needed
    if count_0 == count_1:
        # If counts are equal, we can use either pattern
        return min(mismatches1, mismatches2) // 2
    elif count_0 > count_1:
        # If more '0's, pattern1 must be used
        return mismatches1 // 2
    else:
        # If more '1's, pattern2 must be used
        return mismatches2 // 2

# Given binary string
binary_string = "000101010111101010"
# Calculate minimum swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)