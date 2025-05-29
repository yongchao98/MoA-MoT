def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # Check feasibility based on the length of the string
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps needed for both patterns
    swaps_0_start = 0
    swaps_1_start = 0
    
    for i in range(n):
        if i % 2 == 0:
            if s[i] != '0':
                swaps_0_start += 1
            if s[i] != '1':
                swaps_1_start += 1
        else:
            if s[i] != '1':
                swaps_0_start += 1
            if s[i] != '0':
                swaps_1_start += 1
    
    # Determine the minimum swaps needed
    if n % 2 == 0:
        return min(swaps_0_start // 2, swaps_1_start // 2)
    else:
        if count_0 > count_1:
            return swaps_0_start // 2
        else:
            return swaps_1_start // 2

# Given binary string
binary_string = "100111101010101001"
result = min_swaps_to_alternate(binary_string)
print(result)