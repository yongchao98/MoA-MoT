def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: "010101..."
    mismatches_pattern1 = 0
    # Pattern 2: "101010..."
    mismatches_pattern2 = 0
    
    for i in range(n):
        expected_char_pattern1 = '0' if i % 2 == 0 else '1'
        expected_char_pattern2 = '1' if i % 2 == 0 else '0'
        
        if s[i] != expected_char_pattern1:
            mismatches_pattern1 += 1
        if s[i] != expected_char_pattern2:
            mismatches_pattern2 += 1
    
    # Check if it's possible to make the string alternating
    if mismatches_pattern1 % 2 != 0 and mismatches_pattern2 % 2 != 0:
        return -1
    
    # Calculate the minimum swaps needed
    swaps_pattern1 = mismatches_pattern1 // 2
    swaps_pattern2 = mismatches_pattern2 // 2
    
    # Return the minimum swaps needed
    return min(swaps_pattern1, swaps_pattern2)

# Given binary string
binary_string = "0101010110101100"
result = min_swaps_to_alternate(binary_string)
print(result)