def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: "010101..."
    mismatches1 = 0
    # Pattern 2: "101010..."
    mismatches2 = 0
    
    for i in range(n):
        expected_char1 = '0' if i % 2 == 0 else '1'
        expected_char2 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char1:
            mismatches1 += 1
        if s[i] != expected_char2:
            mismatches2 += 1
    
    # Check if it's possible to make the string alternating
    if mismatches1 % 2 != 0 and mismatches2 % 2 != 0:
        return -1
    
    # Calculate the minimum swaps needed
    swaps1 = mismatches1 // 2
    swaps2 = mismatches2 // 2
    
    if mismatches1 % 2 == 0 and mismatches2 % 2 == 0:
        return min(swaps1, swaps2)
    elif mismatches1 % 2 == 0:
        return swaps1
    else:
        return swaps2

# Given binary string
binary_string = "0101010011110000101111001100"
result = min_swaps_to_alternate(binary_string)
print(result)