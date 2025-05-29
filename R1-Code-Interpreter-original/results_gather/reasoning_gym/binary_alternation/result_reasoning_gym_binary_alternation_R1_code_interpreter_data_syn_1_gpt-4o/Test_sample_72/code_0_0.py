def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between the number of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Pattern 1: "010101..."
    swaps_pattern1 = 0
    for i in range(n):
        expected_char = '0' if i % 2 == 0 else '1'
        if s[i] != expected_char:
            swaps_pattern1 += 1
    
    # Pattern 2: "101010..."
    swaps_pattern2 = 0
    for i in range(n):
        expected_char = '1' if i % 2 == 0 else '0'
        if s[i] != expected_char:
            swaps_pattern2 += 1
    
    # We need to divide by 2 because each swap fixes two positions
    if n % 2 == 0:
        # If length is even, both patterns are possible
        return min(swaps_pattern1 // 2, swaps_pattern2 // 2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return swaps_pattern1 // 2
        else:
            return swaps_pattern2 // 2

# Given binary string
binary_string = "100010011001"
result = min_swaps_to_alternate(binary_string)
print(result)