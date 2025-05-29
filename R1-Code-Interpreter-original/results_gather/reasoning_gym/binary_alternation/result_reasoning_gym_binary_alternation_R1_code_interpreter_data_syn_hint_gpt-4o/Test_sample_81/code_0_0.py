def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of 0s and 1s is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Calculate mismatches for pattern starting with '0' (e.g., "010101...")
    mismatches_0 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            mismatches_0 += 1
    
    # Calculate mismatches for pattern starting with '1' (e.g., "101010...")
    mismatches_1 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            mismatches_1 += 1
    
    # Calculate minimum swaps needed
    if count_0 == count_1:
        # Both patterns are possible
        return min(mismatches_0, mismatches_1) // 2
    elif count_0 > count_1:
        # Only pattern starting with '0' is possible
        return mismatches_0 // 2
    else:
        # Only pattern starting with '1' is possible
        return mismatches_1 // 2

# Given binary string
binary_string = "0001110011"
# Calculate the minimum number of swaps
min_swaps = min_swaps_to_alternate(binary_string)
print(min_swaps)