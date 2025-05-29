def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate swaps for pattern starting with '0' (e.g., "010101...")
    swaps_0 = 0
    swaps_1 = 0
    for i in range(n):
        if i % 2 == 0:
            if s[i] != '0':
                swaps_0 += 1
            if s[i] != '1':
                swaps_1 += 1
        else:
            if s[i] != '1':
                swaps_0 += 1
            if s[i] != '0':
                swaps_1 += 1
    
    # Since each swap fixes two positions, divide by 2
    if count_0 == count_1:
        return min(swaps_0, swaps_1) // 2
    elif count_0 > count_1:
        return swaps_0 // 2
    else:
        return swaps_1 // 2

# Given binary string
binary_string = "100111101010101001"
result = min_swaps_to_alternate(binary_string)
print(result)