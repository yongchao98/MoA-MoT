def min_swaps_to_alternate(s):
    n = len(s)
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # To make the string alternating, we need to swap half of the mismatched characters
    swaps1 = mismatches1 // 2
    swaps2 = mismatches2 // 2
    
    # Check if it's possible to make the string alternating
    if mismatches1 % 2 != 0 and mismatches2 % 2 != 0:
        return -1
    
    return min(swaps1, swaps2)

# Given binary string
binary_string = "0001110110"
result = min_swaps_to_alternate(binary_string)
print(result)