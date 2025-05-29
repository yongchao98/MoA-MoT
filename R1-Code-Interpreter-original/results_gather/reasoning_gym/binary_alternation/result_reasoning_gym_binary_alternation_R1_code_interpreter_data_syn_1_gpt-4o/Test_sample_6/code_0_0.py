def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between count of '0's and '1's is more than 1, it's impossible
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed to match a given pattern
    def calculate_swaps(pattern):
        swaps = 0
        for i in range(n):
            if s[i] != pattern[i % 2]:
                swaps += 1
        return swaps // 2
    
    # Pattern 1: "010101..."
    pattern1 = ['0', '1']
    # Pattern 2: "101010..."
    pattern2 = ['1', '0']
    
    # Calculate swaps for both patterns
    swaps_pattern1 = calculate_swaps(pattern1)
    swaps_pattern2 = calculate_swaps(pattern2)
    
    # If the string length is even, both patterns are possible
    # If the string length is odd, only one pattern is possible
    if n % 2 == 0:
        return min(swaps_pattern1, swaps_pattern2)
    else:
        # If there are more '0's, pattern1 is possible, otherwise pattern2
        if count_0 > count_1:
            return swaps_pattern1
        else:
            return swaps_pattern2

# Given binary string
binary_string = "001101001110100110010101011100"
result = min_swaps_to_alternate(binary_string)
print(result)