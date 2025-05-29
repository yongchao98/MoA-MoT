def min_swaps_to_alternate(s):
    n = len(s)
    
    # Pattern 1: "010101..."
    mismatches_pattern1 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            mismatches_pattern1 += 1
    
    # Pattern 2: "101010..."
    mismatches_pattern2 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            mismatches_pattern2 += 1
    
    # Calculate swaps needed
    if mismatches_pattern1 % 2 != 0 and mismatches_pattern2 % 2 != 0:
        return -1  # Impossible to make alternating
    else:
        swaps_pattern1 = mismatches_pattern1 // 2
        swaps_pattern2 = mismatches_pattern2 // 2
        return min(swaps_pattern1, swaps_pattern2)

# Given binary string
binary_string = "11001001010001111"
result = min_swaps_to_alternate(binary_string)
print(result)