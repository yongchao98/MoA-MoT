def min_swaps_to_alternate(s):
    n = len(s)
    count_0 = s.count('0')
    count_1 = s.count('1')
    
    # If the difference between counts of '0's and '1's is more than 1, return -1
    if abs(count_0 - count_1) > 1:
        return -1
    
    # Function to calculate swaps needed to match a given pattern
    def calculate_swaps(pattern):
        swaps = 0
        for i in range(n):
            if s[i] != pattern[i % 2]:
                swaps += 1
        return swaps // 2
    
    # Possible patterns
    pattern1 = '01' * (n // 2) + '0' * (n % 2)
    pattern2 = '10' * (n // 2) + '1' * (n % 2)
    
    # Calculate swaps for both patterns
    if n % 2 == 0:
        # If length is even, both patterns are possible
        swaps1 = calculate_swaps(pattern1)
        swaps2 = calculate_swaps(pattern2)
        return min(swaps1, swaps2)
    else:
        # If length is odd, only one pattern is possible
        if count_0 > count_1:
            return calculate_swaps(pattern1)
        else:
            return calculate_swaps(pattern2)

# Given binary string
binary_string = "101111001011100000100100101"
result = min_swaps_to_alternate(binary_string)
print(result)