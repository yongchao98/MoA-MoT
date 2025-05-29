def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for both patterns
    mismatches1 = mismatches2 = 0
    for i in range(n):
        if i % 2 == 0:
            if s[i] != '0':
                mismatches1 += 1
            if s[i] != '1':
                mismatches2 += 1
        else:
            if s[i] != '1':
                mismatches1 += 1
            if s[i] != '0':
                mismatches2 += 1
    
    # Calculate minimum swaps needed
    if n % 2 == 0:
        # If length is even, both patterns are possible
        return min(mismatches1 // 2, mismatches2 // 2)
    else:
        # If length is odd, choose the pattern with the majority character
        if count_0 > count_1:
            return mismatches1 // 2
        else:
            return mismatches2 // 2

# Given binary string
binary_string = "010001011100011101001101001010"
# Calculate the minimum number of swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)