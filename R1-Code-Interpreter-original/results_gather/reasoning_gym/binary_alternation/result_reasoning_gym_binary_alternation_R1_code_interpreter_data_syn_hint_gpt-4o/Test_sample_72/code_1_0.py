def min_swaps_to_alternate(s):
    n = len(s)
    
    # Count the number of 0s and 1s
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between the number of 0s and 1s is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
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
    
    # Calculate swaps needed
    swaps_pattern1 = mismatches_pattern1 // 2
    swaps_pattern2 = mismatches_pattern2 // 2
    
    # Determine the minimum swaps needed
    if n % 2 == 0:
        return min(swaps_pattern1, swaps_pattern2)
    else:
        if count_0 > count_1:
            return swaps_pattern1
        else:
            return swaps_pattern2

# Given binary string
binary_string = "100010011001"
# Calculate the minimum swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)