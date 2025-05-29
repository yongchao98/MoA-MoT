def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between the number of 0s and 1s is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Pattern 1: "010101..."
    mismatches_1 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            mismatches_1 += 1
    
    # Pattern 2: "101010..."
    mismatches_2 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            mismatches_2 += 1
    
    # Calculate swaps needed for each pattern
    swaps_1 = mismatches_1 // 2
    swaps_2 = mismatches_2 // 2
    
    # Return the minimum swaps needed
    if count_0 == count_1:
        return min(swaps_1, swaps_2)
    elif count_0 > count_1:
        return swaps_1
    else:
        return swaps_2

# Given binary string
binary_string = "100001101101001011"
result = min_swaps_to_alternate(binary_string)
print(result)