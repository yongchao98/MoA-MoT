def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Generate two possible patterns
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Calculate mismatches for both patterns
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # Calculate minimum swaps
    swaps1 = mismatches1 // 2
    swaps2 = mismatches2 // 2
    
    # Determine the minimum swaps needed
    if count_0 == count_1:
        return min(swaps1, swaps2)
    elif count_0 > count_1:
        return swaps1
    else:
        return swaps2

# Given binary string
binary_string = "1010111011010000"
result = min_swaps_to_alternate(binary_string)
print(result)