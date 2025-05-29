def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed to match a given pattern
    def calculate_swaps(pattern):
        mismatches = 0
        for i in range(n):
            if s[i] != pattern[i]:
                mismatches += 1
        return mismatches // 2
    
    # Generate the two possible patterns
    pattern1 = ''.join(['0' if i % 2 == 0 else '1' for i in range(n)])
    pattern2 = ''.join(['1' if i % 2 == 0 else '0' for i in range(n)])
    
    # Calculate swaps for both patterns
    swaps1 = calculate_swaps(pattern1)
    swaps2 = calculate_swaps(pattern2)
    
    # If the string length is odd, only one pattern is valid
    if n % 2 == 1:
        if count_0 > count_1:
            return swaps1
        else:
            return swaps2
    
    # Return the minimum swaps needed
    return min(swaps1, swaps2)

# Given binary string
binary_string = "00000001111110000111011"
result = min_swaps_to_alternate(binary_string)
print(result)