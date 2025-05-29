def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Pattern 1: Starts with '0'
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    # Pattern 2: Starts with '1'
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Count mismatches for both patterns
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # Calculate minimum swaps needed
    if count_0 == count_1:
        # If counts are equal, we can choose either pattern
        return min(mismatches1 // 2, mismatches2 // 2)
    elif count_0 > count_1:
        # If there are more 0s, we must start with '0'
        return mismatches1 // 2
    else:
        # If there are more 1s, we must start with '1'
        return mismatches2 // 2

# Given binary string
binary_string = "010001011100011101001101001010"
# Calculate the minimum number of swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)