def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between counts of '0's and '1's is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for pattern "010101..."
    swaps_0_start = 0
    swaps_1_start = 0
    
    for i in range(n):
        if i % 2 == 0:  # Even index
            if s[i] != '0':
                swaps_0_start += 1
            if s[i] != '1':
                swaps_1_start += 1
        else:  # Odd index
            if s[i] != '1':
                swaps_0_start += 1
            if s[i] != '0':
                swaps_1_start += 1
    
    # If counts of '0's and '1's are equal, we can choose either pattern
    if count_0 == count_1:
        return min(swaps_0_start, swaps_1_start) // 2
    # If there are more '0's, the pattern must start with '0'
    elif count_0 > count_1:
        return swaps_0_start // 2
    # If there are more '1's, the pattern must start with '1'
    else:
        return swaps_1_start // 2

# Given binary string
binary_string = "010001011100011101001101001010"
result = min_swaps_to_alternate(binary_string)
print(result)