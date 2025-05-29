def min_swaps_to_alternate(s):
    n = len(s)
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    mismatches1 = sum(1 for i in range(n) if s[i] != pattern1[i])
    mismatches2 = sum(1 for i in range(n) if s[i] != pattern2[i])
    
    # If mismatches are odd, it's impossible to swap to that pattern
    if mismatches1 % 2 != 0 and mismatches2 % 2 != 0:
        return -1
    
    # Calculate minimum swaps needed
    swaps1 = mismatches1 // 2 if mismatches1 % 2 == 0 else float('inf')
    swaps2 = mismatches2 // 2 if mismatches2 % 2 == 0 else float('inf')
    
    return min(swaps1, swaps2)

# Given binary string
binary_string = "010110110100"
print(min_swaps_to_alternate(binary_string))